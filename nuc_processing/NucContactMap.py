"""
---- COPYRIGHT ----------------------------------------------------------------

Copyright (C) 20016-2017
Tim Stevens (MRC-LMB) and Wayne Boucher (University of Cambridge)


---- LICENSE ------------------------------------------------------------------

This file is part of NucProcess.

NucProcess is free software: you can redistribute it and/or modify it under the
terms of the GNU Lesser General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

NucProcess is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along
with NucProcess.  If not, see <http://www.gnu.org/licenses/>.


---- CITATION -----------------------------------------------------------------

If you are using this software for academic purposes, we suggest quoting the
following reference:

Stevens TJ, Lando D, Basu S, Atkinson LP, Cao Y, Lee SF, Leeb M, Wohlfahrt KJ,
Boucher W, O'Shaughnessy-Kirwan A, Cramard J, Faure AJ, Ralser M, Blanco E, Morey
L, Sanso M, Palayret MGS, Lehner B, Di Croce L, Wutz A, Hendrich B, Klenerman D,
Laue ED. 3D structures of individual mammalian genomes studied by single-cell
Hi-C. Nature. 2017 Apr 6;544(7648):59-64. doi: 10.1038/nature21429. Epub 2017 Mar
13. PubMed PMID: 28289288; PubMed Central PMCID: PMC5385134.

"""

import numpy as np
import os, sys

from math import ceil
from collections import defaultdict
from .NucSvg import SvgDocument

PROG_NAME = 'nuc_contact_map'
VERSION = '1.1.0'
DESCRIPTION = 'Chromatin contact (NCC format) Hi-C contact map display module for Nuc3D and NucTools'
DEFAULT_BIN_SIZE = 5 # Megabases


def info(msg, prefix='INFO'):
  print('%8s : %s' % (prefix, msg))


def warn(msg, prefix='WARNING'):
  print('%8s : %s' % (prefix, msg))


def fatal(msg, prefix='%s FAILURE' % PROG_NAME):
  print('%8s : %s' % (prefix, msg))
  sys.exit(0)


def _color_func(matrix, colors, zero_color):

  n, m, d = matrix.shape

  color_matrix = np.zeros((n, m, 3), float)
  colors = np.array(colors)
  rgb_0 = np.array(zero_color)

  base_colors = 0.5 * colors

  for k in range(d):
    base_colors[k] += 0.5 * rgb_0

  for i in range(n):
    for j in range(m):
      l = 0.0

      for k in range(d):
        f = matrix[i,j,k]
        if f > 0:
          g = 1.0 - f
          color_matrix[i,j] += (f * colors[k]) + (g * base_colors[k])
          l += 1.0
          break

      if l == 0.0:
        color_matrix[i,j] = rgb_0

      else:
        color_matrix[i,j] /= l

  color_matrix = np.clip(color_matrix, 0, 255)
  color_matrix = np.array(color_matrix, dtype=np.uint8)

  return color_matrix


def _get_trans_dev(trans_counts):

  cp = float(len(trans_counts))

  vals = np.array(list(trans_counts.values()), float)

  if not len(vals):
    return 0.0, '?'

  vals -= vals.min()
  vals /= vals.sum() or 1.0
  vals = vals[vals.argsort()]
  vals = vals.cumsum()

  base = np.arange(0, len(vals))/cp
  deltas = base - vals
  dev = 2.0 * deltas.sum()/cp

  if dev < 0.50:
    score_cat = '>4N'
  elif dev < 0.6:
    score_cat = '4N'
  elif dev < 0.74:
    score_cat = '2N'
  elif dev < 0.76:
    score_cat = '1/2N'
  elif dev < 0.95:
    score_cat = '1N'
  else:
    score_cat = '?'

  return dev, score_cat


def _get_num_isolated(positions, threshold=500000):

  num_isolated = 0
  bin_offsets = ((-1,-1), (-1,0), (-1,1), (0,-1), (0,1), (1,-1), (1,0), (1,1))

  idx = defaultdict(list)

  for i, (p_a, p_b, ag) in enumerate(positions):
    idx[(p_a/threshold, p_b/threshold)].append(i)

  for key in idx:
    if len(idx[key]) == 1:
      pA, pB, ag = positions[idx[key][0]]
      b1, b2 = key

      for j, k in bin_offsets:
        key2 = (b1+j, b2+k)

        if key2 in idx:
          for i2 in idx[key2]:
            pC, pD, ag2 = positions[i2]

            if abs(pC-pA) < threshold and abs(pD-pB) < threshold:
              break

            elif abs(pD-pA) < threshold and abs(pC-pB) < threshold:
              break

          else:
            continue

          break

      else:
        num_isolated += 1

  return num_isolated


def _get_num_isolated_groups(positions, threshold=500000, pos_err=100):

  num_isolated = 0
  found = [0] * len(positions)

  group_dict = defaultdict(list)
  for p_a, p_b, ag in positions:
    group_dict[ag].append((p_a, p_b))

  groups = sorted(group_dict)

  for i, ag_a in enumerate(groups):
    if found[i]: # Already close to something else
      continue

    close = 0
    for j, ag_b in enumerate(groups):
      if j == i:
        continue

      for p_a, p_b in group_dict[ag_a]:
        for p_c, p_d in group_dict[ag_b]:
          if (pos_err < abs(p_c-p_a) < threshold) and (pos_err < abs(p_d-p_b) < threshold):
            close = 1
            found[j] = 1
            break

          elif (pos_err < abs(p_d-p_a) < threshold) and (pos_err < abs(p_c-p_b) < threshold):
            close = 1
            found[j] = 1
            break

        else:
          continue

        break

      if close:
        break

    if not close:
      num_isolated += 1

  return num_isolated


def _get_mito_fraction(contacts, min_sep=1e2, sep_range=(10**6.5, 10**7.5)):

  a, b = sep_range
  in_range = 0
  total = 0

  for chr_pair in contacts:
    chr_a, chr_b = chr_pair

    if chr_a != chr_b:
      continue

    points = np.array(contacts[chr_pair])[:,:2]

    if len(points) < 3:
      continue

    d_seps = np.diff(points, axis=1)
    d_seps = d_seps[(d_seps > min_sep).nonzero()]

    n = len(d_seps)
    smaller = len((d_seps <= a).nonzero()[0])
    larger = len((d_seps >= b).nonzero()[0])

    total += n
    in_range += n-(smaller+larger)

  if not total:
    return 0.0, '?'

  frac = in_range/float(total)

  if frac < 0.30:
    score_cat = 'Non-M'
  elif frac < 0.40:
    score_cat = 'M'
  else:
    score_cat = 'Strong M'

  return frac, score_cat


def load_ncc(ncc_path):

  chromo_limits = {}
  contacts = {}

  with open(ncc_path) as in_file_obj:
    for line in in_file_obj:
      chr_a, f_start_a, f_end_a, start_a, end_a, strand_a, chr_b, f_start_b, f_end_b, start_b, end_b, strand_b, ambig_group, pair_id, swap_pair = line.split()

      f_start_a = int(f_start_a)
      f_end_a = int(f_end_a)
      start_a = int(start_a)
      end_a = int(end_a)
      f_start_b = int(f_start_b)
      f_end_b = int(f_end_b)
      start_b = int(start_b)
      end_b = int(end_b)
      ambig_group = int(ambig_group)

      if chr_a in chromo_limits:
        s, e = chromo_limits[chr_a]
        chromo_limits[chr_a] = (min(s, f_start_a), max(e, f_end_a))
      else:
        chromo_limits[chr_a] = (f_start_a, f_end_a)

      if chr_b in chromo_limits:
        s, e = chromo_limits[chr_b]
        chromo_limits[chr_b] = (min(s, f_start_b), max(e, f_end_b))
      else:
        chromo_limits[chr_b] = (f_start_b, f_end_b)

      # Store data in sorted chromo order
      if chr_a > chr_b:
        chr_a, chr_b = chr_b, chr_a
        start_a, start_b = start_b, start_a
        end_a, end_b = end_b, end_a
        f_start_a, f_start_b = f_start_b, f_start_a
        f_end_a, f_end_b = f_end_b, f_end_a
        strand_a, strand_b = strand_b, strand_a

      chr_pair = (chr_a, chr_b)
      if chr_pair in contacts:
        contact_list = contacts[chr_pair]

      else:
        contact_list = []
        contacts[chr_pair] = contact_list

      if strand_a == '-':
        p_a = f_end_a
      else:
        p_a = f_start_a

      if strand_b == '-':
        p_b = f_end_b
      else:
        p_b = f_start_b

      contact_list.append((p_a, p_b, ambig_group))

  return chromo_limits, contacts


def nuc_contact_map(ncc_path, svg_tag='_contact_map', svg_width=1000, bin_size=5, black_bg=False, color=None,
                    color_ambig=None, font=None, font_size=12, line_width=0.2, min_contig_size=None):

  bin_size = int(bin_size * 1e6)

  if min_contig_size:
    min_contig_size = int(min_contig_size * 1e6)
  else:
    min_contig_size = bin_size

  if svg_tag == '-':
    svg_path = '-'

  else:
    svg_path = '{}{}.svg'.format(os.path.splitext(ncc_path)[0], svg_tag)

  if svg_path and svg_path != '-':
    info('Making contact map for {}'.format(ncc_path))

  chromo_limits, contacts = load_ncc(ncc_path)

  if not chromo_limits:
    fatal('No chromosome contact data read')

  # Get sorted chromosomes, ignore small contigs as appropriate
  chromos = []
  skipped = []
  for chromo in chromo_limits:
    s, e = chromo_limits[chromo]

    if (e-s) < min_contig_size:
      skipped.append(chromo)
      continue

    if chromo.upper().startswith('CHR'):
      c = chromo[3:]
    else:
      c = chromo

    if c.split('.')[-1].upper() in ('A','B'):
      try:
        key = ('%09d' % int(c.split('.')[0]), c.split('.')[-1])
      except ValueError as err:
        key = (c, c.split('.')[-1],)

    else:
      try:
        key = '%09d' % int(c)
      except ValueError as err:
        key = c

    chromos.append((key, chromo))

  if skipped:
    info('Skipped {:,} small chromosomes/contigs < {:,} bp'.format(len(skipped), min_contig_size))

  chromos.sort()
  chromos = [x[1] for x in chromos]

  # Get chromosome matrix index ranges
  grid = []
  chromo_offsets = {}
  chromo_spans = {}
  n = 0
  for chromo in chromos: # In display order
    s, e = chromo_limits[chromo]
    a = int(s/bin_size)
    b = int(ceil(e/float(bin_size))) + 1
    span = b-a
    chromo_offsets[chromo] = s, n # Start bp, start bin index
    chromo_spans[chromo] = span
    n += span
    grid.append(n)# At chromosome edge

  if grid:
    grid.pop() # Don't need last edge

  # Fill contact map matrix, last dim is for (un)ambigous
  data = np.zeros((n, n, 2), float)

  if svg_path and svg_path != '-':
    info('Contact map size %d x %d' % (n, n))

  trans_counts = {}

  n_isol = 0
  n_pairs = 0

  groups = defaultdict(int)
  homolog_groups = set()
  trans_groups = set()
  cis_groups = set()

  for key in contacts:
    for p_a, p_b, ag in contacts[key]:
      groups[ag] += 1

  for i, chr_1 in enumerate(chromos):
    for chr_2 in chromos[i:]:

      if chr_1 > chr_2:
        chr_a, chr_b = chr_2, chr_1
      else:
        chr_a, chr_b = chr_1, chr_2

      contact_list = contacts.get((chr_a, chr_b))

      if contact_list is None: # Nothing for this pair: common for single-cell Hi-C
        if chr_a != chr_b:
          trans_counts[(chr_a, chr_b)] = 0.0

        continue

      ni = _get_num_isolated(contact_list)
      s_a, n_a = chromo_offsets[chr_a]
      s_b, n_b = chromo_offsets[chr_b]

      for p_a, p_b, ag in contact_list:
        if chr_a != chr_b:
          if ('.' in chr_a) and ('.' in chr_b) and (chr_a.split('.')[0] == chr_b.split('.')[0]):
            homolog_groups.add(ag)

          else:
            trans_groups.add(ag)

        else:
          cis_groups.add(ag)

        a = n_a + int((p_a-s_a)/bin_size)
        b = n_b + int((p_b-s_b)/bin_size)

        k = 0 if groups[ag] == 1 else 1

        data[a, b, k] += 1.0
        data[b, a, k] += 1.0

      n_pr = len(contact_list)
      n_isol += ni
      n_pairs += n_pr

      if chr_a != chr_b:
        s1, e1 = chromo_limits[chr_a]
        s2, e2 = chromo_limits[chr_b]
        trans_counts[(chr_a, chr_b)] = (n_pr - ni)/float((e1-s1) * (e2-s2))

  trans_groups -= homolog_groups
  cis_groups -= homolog_groups
  cis_groups -= trans_groups

  n_ambig = len([x for x in groups.values() if x > 1])
  n_homolog = len(homolog_groups)
  n_trans = len(trans_groups)
  n_cis = len(cis_groups)

  n_cont = len(groups)

  ambig_frac = 100.0 * n_ambig / float(n_cont or 1)

  isol_frac = 100.0 * n_isol / float(n_pairs or 1)

  trans_dev, ploidy = _get_trans_dev(trans_counts)

  mito_frac, mito_cat = _get_mito_fraction(contacts)

  stats_text1 = 'Contacts:{:,d} cis:{:,d} trans:{:,d} homolog:{:,d}; ambig:{:.2f}% ; isolated:{:.2f}%'
  stats_text1 = stats_text1.format(n_cont, n_cis, n_trans, n_homolog, ambig_frac, isol_frac)

  stats_text2 = 'ploidy score:{:.2f} ({}) ; mito score:{:.2f} ({})'
  stats_text2 = stats_text2.format(trans_dev, ploidy, mito_frac, mito_cat)

  data = np.log10(data+1.0)

  data[:,:,0] /= data[:,:,0].max() or 1.0
  data[:,:,1] /= data[:,:,1].max() or 1.0

  chromo_labels = []
  for chromo in chromos:
    pos = chromo_offsets[chromo][1] + chromo_spans[chromo]/2

    if chromo.upper().startswith('CHR'):
      chromo = chromo[3:]

    chromo_labels.append((pos, chromo))

  # Make SVG
  offset = 64

  svg_doc = SvgDocument()

  if color:
    color = svg_doc._hex_to_rgb(color)

  if color_ambig:
    color_ambig = svg_doc._hex_to_rgb(color_ambig)

  if black_bg:
    if not color:
      color = [0.0, 200.0, 255.0]
    if not color_ambig:
      color_ambig = [100.0, 0.0, 0.0]

    def color_func(x, c=(color, color_ambig)):
      return _color_func(x, c, [0.0] * 3)

    grid_color = '#303030'

  else:
    if not color:
      color = [0.0, 0.0, 128.0]
    if not color_ambig:
      color_ambig = [180.0, 180.0, 0.0]

    def color_func(x, c=(color, color_ambig)):
      return _color_func(x, c, [255.0] * 3)

    grid_color = '#C0C0C0'

  w = svg_width - 2 * offset
  svg_doc.density_matrix(data, w, w, x_grid=grid, y_grid=grid,
                         x_labels=chromo_labels, y_labels=chromo_labels,
                         x_axis_label='Chromosome', y_axis_label='Chromosome',
                         grid_color=grid_color, font=font,
                         font_size=font_size, line_width=line_width,
                         plot_offset=(offset, offset), color_func=color_func,
                         value_range=None, scale_func=None)

  svg_doc.text(os.path.basename(ncc_path), (offset, 20), size=18)

  svg_doc.text(stats_text1, (offset, 38), size=12)

  svg_doc.text(stats_text2, (offset, 52), size=12)

  if n_ambig:
    x0 = offset + 0.7 * w
    y0 = 52
    x1 = x0 + 8
    y1 = y0 - 8
    svg_doc.rect([x0, y0, x1, y1], color='black', fill='#%02X%02X%02X' % tuple(color))
    svg_doc.text('Unambiguous', (x1+4, y0), size=12)

    x0 = offset + 0.87 * w
    x1 = x0 + 8
    svg_doc.rect([x0, y0, x1, y1], color='black', fill='#%02X%02X%02X' % tuple(color_ambig))
    svg_doc.text('Ambiguous', (x1+4, y0), size=12)

  if svg_path == '-':
    print(svg_doc.svg(svg_width, svg_width))

  elif svg_path:
    svg_doc.write_file(svg_path, svg_width, svg_width)
    info('Saved contact map to %s' % svg_path)

  else:
    return svg_doc.svg(svg_width, svg_width)


def main(argv=None):
  from argparse import ArgumentParser

  if argv is None:
    argv = sys.argv[1:]

  epilog = 'For further help email tjs23@cam.ac.uk or wb104@cam.ac.uk'

  arg_parse = ArgumentParser(prog=PROG_NAME, description=DESCRIPTION,
                             epilog=epilog, prefix_chars='-', add_help=True)

  arg_parse.add_argument(metavar='NCC_FILE', nargs='+', dest='i',
                         help='Input NCC format chromatin contact file(s). Wildcards accepted')

  arg_parse.add_argument('-o', metavar='SVG_FILE_TAG', default='_contact_map',
                         help='Optional name tag to put at end of SVG format contact map file. Use "-" to print SVG to stdout rather than make a file. Default: "_contact_map"')

  arg_parse.add_argument('-w', default=800, metavar='SVG_WIDTH', type=int,
                         help='SVG document width')

  arg_parse.add_argument('-s', default=DEFAULT_BIN_SIZE, metavar='BIN_SIZE', type=int,
                         help='Sequence region size represented by each small square (the resolution) in megabases. Default is %d Mb' % DEFAULT_BIN_SIZE)

  arg_parse.add_argument('-m', default=0.0, metavar='MIN_CONTIG_SIZE', type=float,
                         help='The minimum chromosome/contig sequence length in Megabases for inclusion. Default is the bin size (see -s option)')

  arg_parse.add_argument('-b', default=False, action='store_true',
                         help='Specifies that the contact map should have a black background (default is white)')

  arg_parse.add_argument('-c', nargs=1, metavar='RGB_COLOR',
                         help='Optional main color for the contact points as a 24-bit hexidecimal RBG code e.g. "#0080FF" (with quotes)')

  arg_parse.add_argument('-ca', nargs=1, metavar='RGB_COLOR',
                         help='Optional color for ambigous contact points as a 24-bit hexidecimal RBG code e.g. "#FF0000" (with quotes)')
  args = vars(arg_parse.parse_args(argv))

  ncc_paths = args['i']
  svg_tag = args['o']
  svg_width = args['w']
  bin_size = args['s']
  black_bg = args['b']
  color = args['c']
  color_a = args['ca']
  min_contig_size = args['m']

  if not ncc_paths:
    arg_parse.print_help()
    sys.exit(1)

  if color and isinstance(color, list):
    color = color[0]

  for ncc_path in ncc_paths:
    if not os.path.exists(ncc_path):
      fatal('NCC file could not be found at "{}"'.format(ncc_path))

    nuc_contact_map(ncc_path, svg_tag, svg_width, bin_size, black_bg, color, color_a, min_contig_size=min_contig_size)


if __name__ == "__main__":
  main()
