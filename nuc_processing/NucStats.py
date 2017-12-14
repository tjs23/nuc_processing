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

import json
import os

from .NucProcess import warn, fatal
from .NucContactMap import load_ncc, _get_num_isolated, _get_trans_dev, _get_mito_fraction

PROG_NAME = 'nuc_stats'
VERSION = '1.0.0'
DESCRIPTION = 'Script to aggregate statistics for different single-cell Hi-C datasets from .json and .ncc files'


def nuc_stats(ncc_paths, json_paths, quiet=False, tsv_file_path=None, text_file_path=None):

  paired_files = {}

  for ncc_path in ncc_paths:
    file_root = os.path.splitext(os.path.basename(ncc_path))[0]

    if file_root.endswith('_r_1_2'):
      file_root = file_root[:-6]

    if file_root in paired_files:
      msg = 'NCC file names cannot repeat. File "%s" already considered'
      fatal(msg % paired_files[file_root][0])

    paired_files[file_root] = [ncc_path]

  for json_path in json_paths:
    file_root = os.path.splitext(os.path.basename(json_path))[0]

    if file_root.endswith('_stats'):
      file_root = file_root[:-6]

    if file_root.endswith('_r_1_2'):
      file_root = file_root[:-6]

    if file_root not in paired_files:
      msg = 'JSON file %s cannot be paired with an NCC file. Paired JSON and NCC file names must begin with the same characters.'
      fatal(msg % json_path)

    if len(paired_files[file_root]) == 2:
      msg = 'Cannot pair file "%s". File "%s" already paired with "%s"'
      fatal(msg % (json_path, paired_files[file_root][0], paired_files[file_root][1]))

    paired_files[file_root].append(json_path)

  all_stats = {}

  for file_root in paired_files:
    if len(paired_files[file_root]) == 1:
      msg = 'NCC file %s is not paired with a JSON file and will be ignored.'
      warn(msg % paired_files[file_root][0])
      continue

    ncc_path, json_path = paired_files[file_root]

    with open(json_path) as file_obj:
      stat_dict = json.load(file_obj)

      clip_stats1 = dict(stat_dict['clip_1'])
      clip_stats2 = dict(stat_dict['clip_2'])
      map_1_stats = dict(stat_dict['map_1'])
      map_2_stats = dict(stat_dict['map_2'])

      pair_stats = dict(stat_dict['pair'])
      filter_stats = dict(stat_dict['filter'])
      redundancy_stats = dict(stat_dict['dup'])
      promiscuity_stats = dict(stat_dict['promsic'])
      final_stats = dict(stat_dict['final'])
      frag_sizes = stat_dict['frag_sizes']

      clip_len_1 = clip_stats1['mean_length']
      clip_len_2 = clip_stats2['mean_length']

      n_short_1, n_clip_1 = clip_stats1['too_short']
      n_short_2, n_clip_2 = clip_stats2['too_short']

      n_contacts = final_stats['total_contacts']
      n_cis_near = final_stats['cis_near'][0]
      n_cis_far = final_stats['cis_far'][0]
      n_trans = final_stats['trans'][0]

      n_singe_read, n_ends = redundancy_stats['unique']
      mean_redundancy = redundancy_stats['mean_redundancy']

      n_pairs_vaild = filter_stats['input_pairs'][0]
      n_int_re1 = filter_stats['internal_re1'][0]
      n_adj_re1 = filter_stats['adjacent_re1'][0]
      n_near_cis_pairs, np2 = filter_stats['near_cis_pairs']

      n_other = 0
      for key in ('circular_re1', 'overhang_re1',
                  'too_close', 'too_small', 'too_big',
                  'internal_re2', 'no_end_re2',
                  'low_mappability', 'unknown_contig'):

        n_other += filter_stats[key][0]

      n_unmapped_end, n_pairs = pair_stats['unmapped_end']
      n_unique_map = pair_stats['unique'][0]
      n_anbigous = pair_stats['ambiguous'][0]

      n_promisc, n_prom_inp = promiscuity_stats['promiscuous']
      n_reads_1 = clip_stats1['input_reads']

      n_uniqa, n_readsa = map_1_stats['unique']
      n_uniqb, n_readsb = map_2_stats['unique']
      n_ambiga = map_1_stats['ambiguous'][0]
      n_ambigb = map_2_stats['ambiguous'][0]

      pc_map = (n_uniqa+n_uniqb+n_ambiga+n_ambigb)/(0.01 * (n_readsa+n_readsb))
      pc_short = (n_short_1+n_short_2)/(0.01 * (n_clip_1+n_clip_2))

      n_isol = 0
      n_cont = 1

      chromo_limits, contacts = load_ncc(ncc_path)
      chromos = sorted(chromo_limits)

      trans_counts = {}
      n_cont = 0
      n_isol = 0

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

          nc = len(contact_list)
          n_cont += nc
          n_isol += ni

          if chr_a != chr_b:
            s1, e1 = chromo_limits[chr_a]
            s2, e2 = chromo_limits[chr_b]
            trans_counts[(chr_a, chr_b)] = (nc - ni)/float((e1-s1) * (e2-s2))

      isol_frac = 100.0 * n_isol / float(n_cont or 1)

      trans_dev, ploidy = _get_trans_dev(trans_counts)

      mito_frac, mito_cat = _get_mito_fraction(contacts)

      ccs = '%s/%s' % (ploidy, mito_cat)

      np = 0.01 * float(n_pairs)
      npv = 0.01 * float(n_pairs_vaild)
      ne = 0.01 * float(n_ends)
      npi = 0.01 * float(n_prom_inp)
      nc = 0.01 * float(n_contacts)

      hist, edges = frag_sizes
      nf = sum(hist)
      q1 = 0.25 * nf
      q3 = 0.75 * nf
      r1 = 0
      r2 = 0
      i = 0
      j = 0

      while r1 < q1:
        r1 += hist[i]
        i += 1

      while r2 < q3:
        r2 += hist[j]
        j += 1

      iqr = '%d-%d' % (edges[i], edges[j])
      cl = 0.5*(clip_len_1+clip_len_2)

      all_stats[file_root] = {'rpm':(file_root, cl, n_reads_1, pc_short, pc_map, n_pairs,
                                     n_unique_map/np, n_anbigous/np, n_unmapped_end/np),
                              'qcf':(file_root, n_pairs_vaild, iqr, n_int_re1/npv,
                                     n_adj_re1/npv, n_other/npv, n_singe_read/ne, n_promisc/npi),
                              'pc':(file_root, mean_redundancy, n_contacts, n_cis_near, n_cis_far, n_trans,
                                    n_trans/nc, isol_frac, ccs)}

  cols1 = ('Sample', 'Read length', 'Input pairs', '%Too small', '%Alignment', 'Align pairs', '%Uniq pair', '%Ambig pair')
  cols2 = ('Sample', 'Unique pairs', 'IQ frag size', '%int. RE1', '%adj. RE1', '%other inv', '%one read', '%promisc')
  cols3 = ('Sample', 'Redundancy', 'Total', 'Cis <10 kb', 'Cis >10 kb', 'Trans', '%Trans', '%Iso trans', 'Ploidy/Mito')

  h1 = 'Read Pair Mapping'
  h2 = 'QC & Filtering'
  h3 = 'Processed Contacts'

  if text_file_path or not quiet:
    lines = []
    blank = ' ' * 24

    lines.append(blank + ' | %s' % h1)
    head = ['%-24s' % cols1[0]] + ['%12s' % x for x in cols1[1:]]
    head = ' | '.join(head)

    lines.append(head)
    lines.append('=' * len(head))

    line_fmt = ' | '.join(['{:24s}', '{:12.2f}', '{:12,}', '{:12.2f}', '{:12.2f}', '{:12,}', '{:12.2f}', '{:12.2f}'])

    for file_root in sorted(all_stats):
      line_data = all_stats[file_root]['rpm']

      lines.append(line_fmt.format(*line_data))

    lines.append('-' * len(head))

    lines.append(blank + ' | %s' % h2)
    head = ['%-24s' % cols2[0]] + ['%12s' % x for x in cols2[1:]]
    head = ' | '.join(head)

    lines.append(head)
    lines.append('=' * len(head))

    line_fmt = ' | '.join(['{:24s}', '{:12,}', '{:>12s}', '{:12.2f}', '{:12.2f}', '{:12.2f}', '{:12.2f}', '{:12.2f}'])

    for file_root in sorted(all_stats):
      line_data = all_stats[file_root]['qcf']

      lines.append(line_fmt.format(*line_data))

    lines.append('-' * len(head))

    lines.append(blank + ' | %s' % h3)
    head = ['%-24s' % cols3[0]] + ['%10s' % x for x in cols3[1:]]
    head = ' | '.join(head)

    lines.append(head)
    lines.append('=' * len(head))

    line_fmt = ' | '.join(['{:24s}', '{:10.2f}', '{:10,}', '{:10,}', '{:10,}', '{:10,}', '{:10.2f}', '{:10.2f}', '{:10}'])

    for file_root in sorted(all_stats):
      line_data = all_stats[file_root]['pc']

      lines.append(line_fmt.format(*line_data))

    lines.append('-' * len(head))

    if not quiet:
      for line in lines:
        print(line)

    if text_file_path:
      with open(text_file_path, 'w') as file_obj:
        for line in lines:
          file_obj.write(line + '\n')

  if tsv_file_path:
    with open(tsv_file_path, 'w') as file_obj:
      write = file_obj.write

      write('%s\n' % h1)
      write('\t'.join(cols1) + '\n')

      line_fmt = '\t'.join(['{}', '{:.2f}', '{}', '{:.2f}', '{:.2f}', '{:}', '{:.2f}', '{:.2f}']) + '\n'

      for file_root in sorted(all_stats):
        line_data = all_stats[file_root]['rpm']
        write(line_fmt.format(*line_data))

      write('%s\n' % h2)
      write('\t'.join(cols2) + '\n')

      line_fmt = '\t'.join(['{}', '{:}', '{:s}', '{:.2f}', '{:.2f}', '{:.2f}', '{:.2f}', '{:.2f}']) + '\n'

      for file_root in sorted(all_stats):
        line_data = all_stats[file_root]['qcf']
        write(line_fmt.format(*line_data))

      write('%s\n' % h3)
      write('\t'.join(cols3) + '\n')

      line_fmt = '\t'.join(['{}', '{:.2f}', '{:}', '{:}', '{:}', '{:}', '{:.2f}', '{:.2f}', '{:}']) + '\n'

      for file_root in sorted(all_stats):
        line_data = all_stats[file_root]['pc']
        write(line_fmt.format(*line_data))

# TBD - add totals for whole set?

#  python NucStats.py -n /data/hi-c/GSE80006_zygote/*r_1_2.ncc -j /data/hi-c/GSE80006_zygote/*r_1_2_stats.json -t GSE80006_zygote.txt -o GSE80006_zygote.tsv


if __name__ == '__main__':

  from argparse import ArgumentParser

  epilog = 'Example use: python NucStats.py -n prefix_*.ncc -j prefix_*r_stats.json -t pretty.txt -o stats.tsv'
  epilog += 'For further help email tjs23@cam.ac.uk or wb104@cam.ac.uk'

  arg_parse = ArgumentParser(prog=PROG_NAME, description=DESCRIPTION,
                             epilog=epilog, prefix_chars='-', add_help=True)

  arg_parse.add_argument('-n', nargs='+', metavar='NCC_FILES',
                         help='Input chromosome contact data files in NCC format, as output by nuc_process')

  arg_parse.add_argument('-j', nargs='+', metavar='JSON_FILES',
                         help='Input processing report file in JSON format, as output by nuc_process')

  arg_parse.add_argument('-q', '--quiet', action='store_true', dest='q',
                         help='Supress display of report to screen')

  arg_parse.add_argument('-o', '--out_file',metavar='OUT_TSV_FILE', dest='o',
                         help='Output file path to save report in TSV format')

  arg_parse.add_argument('-t', '--pretty_file',metavar='OUT_TEXT_FILE', dest='t',
                         help='Output file path to save report in a pretty text format')

  args = vars(arg_parse.parse_args())

  json_paths = args['j']
  ncc_paths = args['n']
  quiet = args['q']
  tsv_path = args['o']
  text_path = args['t']

  nuc_stats(ncc_paths, json_paths, quiet, tsv_path, text_path)
