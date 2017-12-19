"""
---- COPYRIGHT ----------------------------------------------------------------

Copyright (C) 20016-2018
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

import sys, os, re, shutil, gzip, json, string
import numpy as np

from collections import defaultdict
from shutil import move
from subprocess import Popen, PIPE, call
from functools import partial
from .NucSvg import SvgDocument
from .NucContactMap import nuc_contact_map
from .common import open_file_r, strip_ext, merge_file_names

PROG_NAME = 'nuc_process'
VERSION = '1.1.1'
DESCRIPTION = 'Chromatin contact paired-read Hi-C processing module for Nuc3D and NucTools'
RE_CONF_FILE = 'enzymes.conf'
RE_SITES = {'MboI'   : '^GATC_',
            'DpnII'  : '^GATC_',
            'AluI'   : 'AG^CT',
            'BglII'  : 'A^GATC_T',
            'HindIII': 'A^AGCT_T'}
QUAL_SCHEMES = ['phred33', 'phred64', 'solexa']
DEFAULT_MIN_QUAL = 10
QUAL_ZERO_ORDS = {'phred33':33, 'phred64':64, 'solexa':64}
FASTQ_READ_CHUNK = 1048576
READ_BUFFER = 2**16
open_file_r = partial(open_file_r, buffering=READ_BUFFER)
MIN_READ_LEN = 20
NUM_MAP_FASTAS = 10
SCORE_TAG = re.compile(r'\sAS:i:(\S+)')
SCORE_TAG_SEARCH = SCORE_TAG.search
NCC_FORMAT = '%s %d %d %d %d %s %s %d %d %d %d %s %d %d %d\n'
LOG_FILE_PATH = None
STAT_FILE_PATH = None
STAT_SECTIONS = {'general':'General',
                 'clip_1':'Clip Reads 1',
                 'clip_2':'Clip Reads 2',
                 'map_1':'Map Reads 1A',
                 'map_2':'Map Reads 2A',
                 'map_3':'Map Reads 1B',
                 'map_4':'Map Reads 2B',
                 'pair':'Read pairing',
                 'filter':'Fragment filtering',
                 'filter_ambig':'Fragment filtering (ambiguous)',
                 'dup':'Duplicate removal',
                 'dup_ambig':'Duplicate removal (ambiguous)',
                 'promsic':'Promiscuous end removal',
                 'promsic_ambig':'Promiscuous end removal (ambiguous)',
                 'final':'Final stats'}

VERBOSE = True
SESSION_KEY = ''.join(np.random.choice(tuple(string.ascii_letters), 8))
TEMP_EXT = '_temp_%s' % (SESSION_KEY)
INTERRUPTED = False

if os.path.exists(RE_CONF_FILE):
  for line in open(RE_CONF_FILE):
    name, site = line.split()
    RE_SITES[name] = site


def compress_file(file_path):

  in_file_obj = open(file_path, 'rU', READ_BUFFER)
  out_file_path = file_path + '.gz'

  out_file_obj = gzip.open(out_file_path, 'wb')
  out_file_obj.writelines(in_file_obj)

  in_file_obj.close()
  out_file_obj.close()

  os.unlink(file_path)

  return out_file_path


def tag_file_name(file_path, tag, file_ext=None, sep='_', ncc_tag='_nuc'):

  if tag:
    tag = sep + tag

  dir_name, file_name = os.path.split(file_path)

  if file_name.endswith('.gz'):
    file_root, file_ext_old = os.path.splitext(file_name[:-3])
    file_name = file_root + tag + (file_ext or file_ext_old) + '.gz'

  elif file_name.endswith('.ncc') or (file_ext == '.ncc'):
    file_root, file_ext_old = os.path.splitext(file_name)

    if '_ambig' in file_root:
      ambig_tag = '_ambig'
    else:
      ambig_tag = ''

    if ncc_tag in file_root:
      file_root = file_root[:file_root.rindex(ncc_tag)]

    file_name = file_root + ncc_tag + tag + ambig_tag + (file_ext or file_ext_old)

  else:
    file_root, file_ext_old = os.path.splitext(file_name)
    file_name = file_root + tag + (file_ext or file_ext_old)

  file_path = os.path.join(dir_name, file_name)

  return file_path


def write_sam_file(ncc_file_path, ref_sam_file_1, ref_sam_file_2):
   """
   Convert the NCC text format to a SAM file by refering back to the original mapped SAM files.
   Assumes entries in original SAM files are paired, for the moment
   """

   sam_file_path = os.path.splitext(ncc_file_path)[0] + '.sam'
   ncc_file_obj = open_file_r(ncc_file_path)
   sam_file_obj = open(sam_file_path, 'w')

   sam_ids = {}
   for line in ncc_file_obj:
     pair_id, swap_pair = line.split()[-2:]
     sam_ids[int(pair_id)] = int(swap_pair)

   in_file_obj_1 = open_file_r(ref_sam_file_1)
   in_file_obj_2 = open_file_r(ref_sam_file_2)

   read2 = in_file_obj_2.readline
   write = sam_file_obj.write

   for line1 in in_file_obj_1:
     if line1[0] == '@':
       write(line1)
     else:
       write('@PG\tNucProcess {}\n'.format(VERSION))
       break

   in_file_obj_1.seek(0)

   # Skip to end of header
   line2 = read2()
   while line2:
     if line2[0] != '@':
       break

     line2 = read2()

   for line1 in in_file_obj_1:
     if line1[0] == '@':
       continue

     id1 = int(line1[:10])
     id2 = id1-1

     while (id2 < id1) and line2:
       id2 = int(line2[:10])
       line2 = read2()

     if (id1 == id2) and (id1 in sam_ids):
       if sam_ids[id1]:
         line2, line1 = line1, line2

       write('%s\n%s\n' % pair_sam_lines(line1[10:], line2[10:]))

   return sam_file_path


def remove_promiscuous(ncc_file, num_copies=1, keep_files=True, zip_files=False,
                       resolve_limit=1e3, close_cis=1e4, ambig=False):
  """
  Function to remove contacts with promiscuous ends from an NCC format file.
  Promiscuous ends occur where a specififc restriction fragment evd is involved
  in more contacts that would normally be allowed given the ploidy of the cell.

  resolve_limit : Allow two suitably close promiscous ends if the pairs are long range cis or trans
  """

  clean_ncc_file = tag_file_name(ncc_file, 'clean')
  clean_ncc_file_temp = clean_ncc_file + TEMP_EXT

  if INTERRUPTED and os.path.exists(clean_ncc_file) and not os.path.exists(clean_ncc_file_temp):
    return clean_ncc_file

  if keep_files:
    promiscous_ncc_file = tag_file_name(ncc_file, 'promisc')
    promisc_file_obj = open(promiscous_ncc_file, 'w')
    write_promisc = promisc_file_obj.write

  in_file_obj = open(ncc_file)
  clean_file_obj = open(clean_ncc_file_temp, 'w')
  write_clean = clean_file_obj.write

  frag_counts = defaultdict(set)
  n_promiscuous = 0
  n_resolved = 0
  n_clean = 0

  for line in in_file_obj:
    line_data = line.split()
    chr_a = line_data[0]
    chr_b = line_data[6]
    f_start_a = int(line_data[1])
    f_start_b = int(line_data[7])
    strand_a = line_data[5]
    strand_b = line_data[11]
    ambig_group = int(line_data[12])
    frag_counts[(chr_a, f_start_a, strand_a)].add((chr_b, f_start_b, strand_b, ambig_group))
    frag_counts[(chr_b, f_start_b, strand_b)].add((chr_a, f_start_a, strand_a, ambig_group))

  remove = set()
  remove_ambig_group = set()

  if hasattr(frag_counts, 'iteritems'):
    items = frag_counts.iteritems
  else:
    items = frag_counts.items

  for re_start, re_ends in items():
    if len(re_ends) > num_copies:
      ambig_groups = {x[3] for x in re_ends}

      if len(ambig_groups) == 1: # Multiple ends are due to mapping ambiguity only
        continue

      if len(re_ends) == 2:
        chr_a, f_start_a, strand_a = re_start
        chr_b1, f_start_b1, strand_b1, ambig_b1 = re_ends.pop()
        chr_b2, f_start_b2, strand_b2, ambig_b2 = re_ends.pop()

        if (chr_b1 == chr_b2) and abs(f_start_b1-f_start_b2) < resolve_limit: # Other ends are very close to each other
          if chr_a != chr_b1: # Trans to this end
            n_resolved += 1

          elif abs(f_start_a-f_start_b1) > close_cis and abs(f_start_a-f_start_b2) > close_cis: # Far from this end
            n_resolved += 1

          else: # Other end too close to this end
            remove.add(re_start)
            remove_ambig_group.add(ambig_b1)
            remove_ambig_group.add(ambig_b2)

        else: # Other ends too far apart
          remove.add(re_start)
          remove_ambig_group.add(ambig_b1)
          remove_ambig_group.add(ambig_b2)

      else: # Cannot resolve with more than two ends
        remove.add(re_start)
        for re_end in re_ends: # All associated ambiguity groups annuled
          remove_ambig_group.add(re_end[3])

  in_file_obj.seek(0)
  for line in in_file_obj:
    line_data = line.split()
    chr_a = line_data[0]
    chr_b = line_data[6]
    f_start_a = int(line_data[1])
    f_start_b = int(line_data[7])
    strand_a = line_data[5]
    strand_b = line_data[11]
    ambig_group = line_data[12]

    if ((ambig_group in remove_ambig_group)
        or ((chr_a, f_start_a, strand_a) in remove)
        or ((chr_b, f_start_b, strand_b) in remove)):
      n_promiscuous += 1

      if keep_files:
        write_promisc(line)

    else:
      n_clean += 1
      write_clean(line)

  in_file_obj.close()

  if keep_files:
    if zip_files:
      compress_file(ncc_file)
      compress_file(promiscous_ncc_file)

  else:
    os.unlink(ncc_file)

  n = n_promiscuous + n_clean
  n_clean -= n_resolved

  stats = [('input_pairs', n),
           ('clean',(n_clean, n)),
           ('promiscuous',(n_promiscuous, n)),
           ('resolved',(n_resolved, n)),
           ('accepted',(n_clean+n_resolved, n)),
           ]

  stat_key = 'promsic_ambig' if ambig else 'promsic'
  log_report(stat_key, stats)

  move(clean_ncc_file_temp, clean_ncc_file)

  return clean_ncc_file


def get_ncc_stats(ncc_file, hom_chromo_dict, far_min=100000):

  n = 0
  n_trans = 0
  n_cis_near = 0
  n_cis_far = 0

  groups = defaultdict(int)
  homolog_groups = set()
  trans_groups = set()
  near_groups = set()
  far_groups = set()

  with open(ncc_file) as in_file_obj:
    for line in in_file_obj:
      line_data = line.split()
      chr_a = line_data[0]
      chr_b = line_data[6]
      group = int(line_data[12])
      groups[group] += 1
      n += 1

      if chr_a != chr_b:
        if chr_b == hom_chromo_dict.get(chr_a):
          homolog_groups.add(group)
        else:
          trans_groups.add(group)

      else:
        strand_a = line_data[5]
        strand_b = line_data[11]

        if strand_a == '+':
          p1 = int(line_data[2]) # Ligation junction is ahead at 3' end of RE frag
        else:
          p1 = int(line_data[1])

        if strand_b == '+':
          p2 = int(line_data[8])
        else:
          p2 = int(line_data[7])

        if abs(p1-p2) < far_min:
          near_groups.add(group)
        else:
          far_groups.add(group)

  n_ambig = len([x for x in groups.values() if x > 1])

  trans_groups -= homolog_groups

  near_groups -= homolog_groups
  near_groups -= trans_groups

  far_groups -= homolog_groups
  far_groups -= trans_groups
  far_groups -= near_groups

  n_homolog = len(homolog_groups)
  n_trans = len(trans_groups)
  n_cis_near = len(near_groups)
  n_cis_far = len(far_groups)
  n_contacts = n_homolog + n_trans + n_cis_near + n_cis_far

  stats = [('total_pairs', n),
           ('total_contacts', n_contacts),
           ('ambiguous_contacts', (n_ambig, n_contacts)),
           ('cis_homolog',(n_homolog, n_contacts)),
           ('cis_near',(n_cis_near, n_contacts)),
           ('cis_far',(n_cis_far, n_contacts)),
           ('trans',(n_trans, n_contacts))]

  return stats


def remove_redundancy(ncc_file, min_repeats=2, keep_files=True, zip_files=False, use_re_fragments=True, ambig=False):

  """
  A:B B:A redundancey taken care of at this stage because the NCC file pairs are internally sorted

  Option to remove redundancy at he RE fragment level (or otherwise at the read level)

  Choose the longest read (not most common?) as the representitive for a merge
  """

  out_file_name = tag_file_name(ncc_file, 'multi_read')
  out_file_name_temp = out_file_name + TEMP_EXT

  if INTERRUPTED and os.path.exists(out_file_name) and not os.path.exists(out_file_name_temp):
    return out_file_name

  # First make temporary sorted file
  sort_file_name = tag_file_name(ncc_file, 'sort')
  cmd_args = ['sort', ncc_file]
  call(cmd_args, shell=False, stderr=None, stdin=None, stdout=open(sort_file_name, 'w'))

  sort_file_obj = open(sort_file_name, 'r')
  out_file_obj = open(out_file_name_temp, 'w')

  n_unique = 0
  n_redundant = 0
  n_pairs = 0
  mean_redundancy = 0.0

  # Compare using strand information given we want to know which fragment END is used
  if use_re_fragments:
    ncc_idx = (0,1,5,6,7,11) # Compare on RE fragment: chr_a, f_start_a, strand_a, chr_b, f_start_b, strand_b
  else:
    ncc_idx = (0,3,5,6,9,11) # Compare on read starts - does not consider ends as these can vary in read replicates

  # Calculate sizes of ambiguity groups

  nc = 0
  ambig_sizes = defaultdict(int)
  for line in sort_file_obj:
    ambig_group = line.split()[12]
    ambig_sizes[ambig_group] += 1
    nc += 1

  # Remove repeats

  excluded_groups = set()
  group_reps = defaultdict(int)
  sort_file_obj.seek(0)
  line_keep = sort_file_obj.readline()
  write = out_file_obj.write

  if line_keep: # Could be empty
    line_data = line_keep.split()

    n = 1
    keep_data = [line_data[i] for i in ncc_idx]
    len_keep = abs(int(line_data[4]) - int(line_data[3])) + abs(int(line_data[10]) - int(line_data[9]))
    ambig_keep = line_data[12]
    n_ambig_keep = ambig_sizes[ambig_keep]

    # Remove redundancy, keeping the least ambigous or else the longest
    for line in sort_file_obj:
      n_pairs += 1
      line_data = line.split()
      curr_data = [line_data[i] for i in ncc_idx]
      len_curr = abs(int(line_data[4]) - int(line_data[3])) + abs(int(line_data[10]) - int(line_data[9]))
      ambig_curr = line_data[12]
      n_ambig_curr = ambig_sizes[ambig_curr]

      if curr_data == keep_data:
        n += 1

        if n_ambig_curr == n_ambig_keep: # Equally ambiguous
          if len_curr > len_keep: # This, longer repeat is better
            excluded_groups.add(ambig_keep) # Remove prev kept group
            len_keep = len_curr
            line_keep = line
            n_ambig_keep = n_ambig_curr
            ambig_keep = ambig_curr

          else: # This repeat is worse
            excluded_groups.add(ambig_curr) # Remove this group

        elif n_ambig_curr < n_ambig_keep: # Better to keep this less ambiguous repeat
          excluded_groups.add(ambig_keep) # Remove prev kept group
          len_keep = len_curr
          line_keep = line
          n_ambig_keep = n_ambig_curr
          ambig_keep = ambig_curr

        else: # Keep the old one
          excluded_groups.add(ambig_curr) # Remove this group

        continue

      else: # Not a repeat, write previous
        group_reps[ambig_keep] += n
        write(line_keep)
        mean_redundancy += n

        line_keep = line
        len_keep = len_curr
        keep_data = curr_data
        ambig_keep = ambig_curr
        n_ambig_keep = n_ambig_curr
        n = 1

    # Write remaining
    write(line_keep)

  out_file_obj.close()

  if keep_files:
    if zip_files:
      compress_file(ncc_file)

  else:
    os.unlink(ncc_file)

  os.unlink(sort_file_name) # Remove temp sorted file

  # Remove excluded ambig groups and non-supported

  if keep_files:
    uniq_file_name = tag_file_name(ncc_file, 'unique_read')
    uniq_file_obj = open(uniq_file_name, 'w')
    uniq_write = uniq_file_obj.write

  n_excluded = 0
  out_file_obj = open(out_file_name, 'w')
  with open(out_file_name_temp, 'r') as temp_file_obj:
    write = out_file_obj.write

    for line in temp_file_obj:
      ambig = line.split()[12]

      if ambig in excluded_groups:
        n_excluded += 1

      elif group_reps[ambig]/float(ambig_sizes[ambig]) < min_repeats:
        if keep_files:
          uniq_write(line)
        n_unique += 1

      else:
        write(line)
        n_redundant += 1

  out_file_obj.close()

  if keep_files:
    uniq_file_obj.close()

  n = n_unique + n_redundant
  mean_redundancy /= float(n) or 1.0

  stats = [('input_pairs', n_pairs),
           ('defunct_ambig_group', (n_excluded, n)),
           ('unique', (n_unique, n)),
           ('redundant', (n_redundant, n)),
           ('mean_redundancy', mean_redundancy)]

  stat_key = 'dup_ambig' if ambig else 'dup'
  log_report(stat_key, stats)

  os.unlink(out_file_name_temp)

  return out_file_name


def filter_pairs(pair_ncc_file, re1_files, re2_files, sizes=(100,2000), keep_files=True, zip_files=False,
                 chromo_name_dict=None, min_mappability=1, re2_tolerance=1000, too_close_limit=1000, ambig=False):

  filter_file = tag_file_name(pair_ncc_file, 'filter_accepted', '.ncc')
  filter_file_temp = filter_file + TEMP_EXT

  if INTERRUPTED and os.path.exists(filter_file) and not os.path.exists(filter_file_temp):
    return filter_file, {}

  re1_frag_dict = {}
  re1_end_dict = {}
  chromo_name_dict1 = {}
  re2_frag_dict = {}
  re2_end_dict = {}
  chromo_name_dict2 = {}

  for re1_file in re1_files:
    frags, ends, names = read_re_frag_file(re1_file)
    re1_frag_dict.update(frags)
    re1_end_dict.update(ends)
    chromo_name_dict1.update(names)

  for re2_file in re2_files:
    frags, ends, names = read_re_frag_file(re2_file)
    re2_frag_dict.update(frags)
    re2_end_dict.update(ends)
    chromo_name_dict2.update(names)

  if not chromo_name_dict:
    chromo_name_dict = {}

  out_file_objs = {}
  write_funcs = {}
  out_file_names = {}
  counts = {}
  size_counts = defaultdict(int)
  min_size, max_size = sorted(sizes)
  frag_sizes = []

  for tag in ('accepted', 'too_small', 'too_big', 'circular_re1', 'internal_re1',
              'internal_re2', 'no_end_re2', 'overhang_re1', 'too_close',
              'adjacent_re1', 'unknown_contig', 'excluded_ambig_group'):
    counts[tag] = 0

    if keep_files or (tag == 'accepted'):
      if tag == 'accepted':
        file_name = filter_file_temp
        file_obj = open(filter_file_temp, 'w')
      else:
        file_name = tag_file_name(pair_ncc_file, 'filter_' + tag, '.ncc')
        file_obj = open(file_name, 'w')

      out_file_objs[tag] = file_obj
      out_file_names[tag] = file_name
      write_funcs[tag] = file_obj.write

  counts['near_cis_pairs'] = 0
  counts['far_cis_pairs'] = 0
  counts['trans_pairs'] = 0
  counts['input_pairs'] = 0
  counts['excluded_group_member'] = 0

  def count_write(tag, line):
    counts[tag] += 1

    if keep_files or (tag == 'accepted'):
      write_funcs[tag](line)

  in_file_obj = open_file_r(pair_ncc_file)
  excluded_groups = set()

  for line in in_file_obj:
    # NCC file starts are always the start of the READ, not what was in the BAM file
    # which means start > end where strand is "-"

    chr_a, f_start_a, f_end_a, start_a, end_a, strand_a, chr_b, f_start_b, f_end_b, start_b, end_b, strand_b, ambig_group, pair_id, swap_pair = line.split()

    counts['input_pairs'] += 1

    pos_strand_a = strand_a == '+'
    pos_strand_b = strand_b == '+'

    start_a = int(start_a)
    start_b = int(start_b)

    end_a = int(end_a)
    end_b = int(end_b)

    ambig_group = int(ambig_group)

    # Convert contig code names to proper chromosome names where possible
    chr_a = chromo_name_dict1.get(chr_a, chr_a)
    chr_b = chromo_name_dict1.get(chr_b, chr_b)

    chr_name_a = chromo_name_dict.get(chr_a, chr_a)
    chr_name_b = chromo_name_dict.get(chr_b, chr_b)

    if chr_a not in re1_end_dict:
      count_write('unknown_contig', line)
      continue

    if chr_b not in re1_end_dict:
      count_write('unknown_contig', line)
      continue

    # Find which RE fragments the read positions wre within
    # index of frag with pos immediately less than end
    # Note: modulo is due to circular chromosomes
    re1_a_idx = np.searchsorted(re1_end_dict[chr_a], [start_a])[0] % len(re1_end_dict[chr_a])
    re1_b_idx = np.searchsorted(re1_end_dict[chr_b], [start_b])[0] % len(re1_end_dict[chr_b])

    re1_a_start, re1_a_end, mappability_a = re1_frag_dict[chr_a][re1_a_idx]
    re1_b_start, re1_b_end, mappability_b = re1_frag_dict[chr_b][re1_b_idx]

    # record fragment lengths...

    if pos_strand_a:
      delta_re1_a = re1_a_end - start_a # separation from ligation junction
      p1, p2 = start_a, end_a # sorted GENOME positions
    else:
      delta_re1_a = start_a - re1_a_start
      p1, p2 = end_a, start_a

    if pos_strand_b:
      delta_re1_b = re1_b_end - start_b
      p3, p4 = start_b, end_b
    else:
      delta_re1_b = start_b - re1_b_start
      p3, p4 = end_b, start_b

    size_t = delta_re1_a + delta_re1_b

    size_counts[size_t] += 1
    frag_sizes.append(size_t)

    size_ok = min_size <= size_t <= max_size

    # Add RE fragment positions
    line = NCC_FORMAT % (chr_name_a, start_a, end_a, re1_a_start, re1_a_end, strand_a,
                         chr_name_b, start_b, end_b, re1_b_start, re1_b_end, strand_b,
                         ambig_group, int(pair_id), int(swap_pair))

    if re2_files:
      # Check Read is at the end of different RE2 fragments
      # Note: modulo is due to circular chromosomes
      re2_a_idx = np.searchsorted(re2_end_dict[chr_a], [start_a])[0] % len(re2_end_dict[chr_a])
      re2_b_idx = np.searchsorted(re2_end_dict[chr_b], [start_b])[0] % len(re2_end_dict[chr_b])

      # Internal RE2 not necessarily a problem if reads are in different RE1 fragements, e.g. potential ligation junction within
      if re2_a_idx == re2_b_idx and chr_a == chr_b and re1_a_idx == re1_b_idx:
        count_write('internal_re2', line)
        excluded_groups.add(ambig_group) # The real pair could definitely be a duff one
        continue

      re2_a_start, re2_a_end, mappability_re2_a = re2_frag_dict[chr_a][re2_a_idx]
      re2_b_start, re2_b_end, mappability_re2_b = re2_frag_dict[chr_b][re2_b_idx]

      if pos_strand_a:
        delta_re2_a = start_a - re2_a_start  # separation from 5'
      else:
        delta_re2_a = re2_a_end - start_a

      if pos_strand_b:
        delta_re2_b = start_b - re2_b_start
      else:
        delta_re2_b = re2_b_end - start_b

      if abs(delta_re2_a) > re2_tolerance:
        count_write('no_end_re2', line)
        excluded_groups.add(ambig_group) # The real pair could definitely be a duff one
        continue

      if abs(delta_re2_b) > re2_tolerance:
        count_write('no_end_re2', line)
        excluded_groups.add(ambig_group) # The real pair could definitely be a duff one
        continue

    if chr_a == chr_b:

      if re1_a_idx == re1_b_idx: # Same Re1 fragment
        if start_a < start_b:
          if pos_strand_b and not pos_strand_a: # Reads go outwards in same frag: <-A B->
            count_write('circular_re1', line)
            excluded_groups.add(ambig_group)
            continue

        else:
          if pos_strand_a and not pos_strand_b: # Reads go outwards in same frag: <-B A->
            count_write('circular_re1', line)
            excluded_groups.add(ambig_group)
            continue

        if (p1 < re1_a_start < p2) or (p1 < re1_a_end < p2): # Read A overshoots Re1 fragment end
          count_write('overhang_re1', line)

        elif (p3 < re1_b_start < p4) or (p3 < re1_b_end < p4): # Read B overshoots Re1 fragment end
          count_write('overhang_re1', line)

        else:
          count_write('internal_re1', line)

        excluded_groups.add(ambig_group) # The real pair could definitely be a duff one

      elif abs(re1_a_idx-re1_b_idx) < 2:
        count_write('adjacent_re1', line) # Mostly re-ligation
        excluded_groups.add(ambig_group) # The real pair could definitely be useless

      else: # Different re1 fragment

        if pos_strand_a != pos_strand_b:

          if pos_strand_a and (p1 < p3): # Sequencing toward each other
            delta = p4 - p1 # separation of pair

            if delta < too_close_limit:
              count_write('too_close', line)  # Pair is possible even without ligation
              excluded_groups.add(ambig_group)  # The real pair could definitely be useless
              continue

          elif pos_strand_b and (p3 < p1): # Sequencing toward each other
            delta = p2 - p3

            if delta < too_close_limit:
              count_write('too_close', line)  # Pair is possible even without ligation
              excluded_groups.add(ambig_group) # The real pair could definitely be useless
              continue

        if size_ok:
          if min(abs(p1-p3), abs(p2-p4)) < 1e4:
            counts['near_cis_pairs'] += 1
          else:
            counts['far_cis_pairs'] += 1

          count_write('accepted', line)

        elif size_t < min_size:
          count_write('too_small', line)
          excluded_groups.add(ambig_group)

        else:
          count_write('too_big', line)
          excluded_groups.add(ambig_group)

    else: # Good stuff, trans contact, many errors impossible
      if size_ok:
        counts['trans_pairs'] += 1
        count_write('accepted', line)

      elif size_t < min_size:
        count_write('too_small', line)
        excluded_groups.add(ambig_group)

      else:
        count_write('too_big', line)
        excluded_groups.add(ambig_group)

  # Remove complete excluded ambiguity groups at the end:
  # - Removes a whole ambiguity group if only one possibilty was suspect as it could be the real contact

  out_file_objs['accepted'].close()
  del out_file_objs['accepted']

  out_file_obj = open(filter_file, 'w')
  with open_file_r(out_file_names['accepted']) as file_obj:
    write = out_file_obj.write

    for line in file_obj:
      ambig_group = int(line.split()[12])

      if ambig_group in excluded_groups:
        count_write('excluded_ambig_group', line)

      else:
        write(line)

  counts['accepted'] -= counts['excluded_ambig_group']

  for tag in out_file_objs:
    out_file_objs[tag].close()

  n = counts['input_pairs']
  stats_list = []
  for key in ('input_pairs', 'accepted', 'near_cis_pairs', 'far_cis_pairs', 'trans_pairs',
              'internal_re1', 'adjacent_re1', 'circular_re1', 'overhang_re1',
              'too_close', 'too_small', 'too_big', 'excluded_ambig_group',
              'internal_re2', 'no_end_re2',
              'unknown_contig'):

    stats_list.append((key, (counts[key], n))) # So percentages can be calculated in reports

  del out_file_names['accepted'] # Main output never compressed at this stage

  if keep_files and zip_files:
    for tag in out_file_objs:
      compress_file(out_file_names[tag])

  frag_hist, size_edges = np.histogram(frag_sizes, bins=50, range=(0, sizes[1]))

  # JSON encoder in Python 3 has become picky about NumPy datatypes
  frag_hist = [int(x) for x in frag_hist]
  size_edges = [float(x) for x in size_edges]

  stat_key = 'filter_ambig' if ambig else 'filter'
  frag_key = 'frag_sizes_ambig' if ambig else 'frag_sizes'

  log_report(stat_key, stats_list, {frag_key:(frag_hist, size_edges)})

  os.unlink(filter_file_temp)

  return filter_file, out_file_names


def pair_sam_lines(line1, line2):
  """
  Convert from single to paired SAM format
  """

  read1 = line1.split('\t')
  read2 = line2.split('\t')

  # Relevant bitwise flags (flag in an 11-bit binary number)
  # 1 The read is one of a pair
  # 2 The alignment is one end of a proper paired-end alignment
  # 4 The read has no reported alignments
  # 8 The read is one of a pair and has no reported alignments
  # 16 The alignment is to the reverse reference strand
  # 32 The other mate in the paired-end alignment is aligned to the reverse reference strand
  # 64 The read is the first (#1) mate in a pair
  # 128 The read is the second (#2) mate in a pair

  # The reads were mapped as single-end data, so should expect flags of
  # 0 (map to the '+' strand) or 16 (map to the '-' strand)
  # Output example: a paired-end read that aligns to the reverse strand
  # and is the first mate in the pair will have flag 83 (= 64 + 16 + 2 + 1)

  bits1 = int(read1[1])
  bits2 = int(read2[1])

  # Ignore non-alignments (and reads rejected due to -m 1 in Bowtie)
  if (bits1 & 0x4) or (bits2 & 0x4):
    return (0, 0)

  # The flag should now indicate this is paired-end data
  bits1 |= 0x1
  bits1 |= 0x2
  bits2 |= 0x1
  bits2 |= 0x2

  # Indicate if the pair is on the reverse strand
  if bits1 & 0x10:
    bits2 |= 0x20

  if bits2 & 0x10:
    bits1 |= 0x20

  # Determine if this is the first or the second pair
  bits1 |= 0x40
  bits2 |= 0x80

  # Insert the modified bitwise flags into the reads
  read1[1] = bits1
  read2[1] = bits2

  # Determine the RNEXT and PNEXT values (i.e. the positional values of a read's pair)
  # RNEXT
  if read1[2] == read2[2]:
    read1[6] = '='
    read2[6] = '='
  else:
    read1[6] = read2[2]
    read2[6] = read1[2]

  # PNEXT
  read1[7] = read2[3]
  read2[7] = read1[3]

  line1 = '\t'.join([str(t) for t in read1])
  line2 = '\t'.join([str(t) for t in read2])

  return line1, line2


def _sort_sam_file_ids(sam_file_path, buf_size=2**16):
  # This function only needed if Bowtie2 does not use --reorder option but does use multiple cores (-p option)

  sort_sam_file_name = tag_file_name(sam_file_path, 'sort', '.sam')

  temp_main = tag_file_name(sam_file_path, 'temp1', '.sam')
  temp_sort = tag_file_name(sam_file_path, 'temp2', '.sam')

  in_file_obj = open(sam_file_path, 'r', buf_size)
  out_file_obj_h = open(sort_sam_file_name, 'w', buf_size)
  out_file_obj_m = open(temp_main, 'w', buf_size)
  write = out_file_obj_m.write

  for line in in_file_obj:
    if line[0] == '@':
      out_file_obj_h.write(line)
    else:
      write(line)

  in_file_obj.close()
  out_file_obj_h.close()
  out_file_obj_m.close()

  call(['sort', temp_main], shell=False, stderr=None, stdin=None, stdout=open(temp_sort, 'w'))

  in_file_obj = open(temp_sort)
  out_file_obj = open(sort_sam_file_name, 'a')
  out_file_obj.writelines(in_file_obj)
  out_file_obj.close()
  in_file_obj.close()

  os.unlink(temp_main)
  os.unlink(temp_sort)

  return sort_sam_file_name


def _skip_sam_header(readline):

  line = readline()

  while line and line[0] == '@':
    line = readline()

  return line


def _write_ncc_line(ncc_data_a, ncc_data_b, ambig_group, _id, write_func):
  # write NCC line in chromosome, position order for merging etc

  chr_a = ncc_data_a[0]
  chr_b = ncc_data_b[0]
  start_a = ncc_data_a[3]
  start_b = ncc_data_b[3]

  if chr_a == chr_b:
    if start_a > start_b:
      line_data = ncc_data_b + ncc_data_a + (ambig_group, _id, 1)

    else:
      line_data = ncc_data_a + ncc_data_b + (ambig_group, _id, 0)

  elif chr_a > chr_b:
    line_data = ncc_data_b + ncc_data_a + (ambig_group, _id, 1)

  else:
    line_data = ncc_data_a + ncc_data_b + (ambig_group, _id, 0)

  write_func(NCC_FORMAT % line_data)


def pair_mapped_hybrid_seqs(sam_file1, sam_file2, sam_file3, sam_file4, unpair_path1, unpair_path2, file_root, ambig=True, unique_map=False, max_unpairable=1):

  paired_ncc_file_name = tag_file_name(file_root, 'pair', '.ncc')
  ambig_ncc_file_name = tag_file_name(file_root, 'pair_ambig', '.ncc')
  paired_ncc_file_name_temp = paired_ncc_file_name + TEMP_EXT
  ambig_ncc_file_name_temp = ambig_ncc_file_name + TEMP_EXT

  if (INTERRUPTED and os.path.exists(paired_ncc_file_name)
      and os.path.exists(ambig_ncc_file_name)
      and not os.path.exists(paired_ncc_file_name_temp)
      and not os.path.exists(ambig_ncc_file_name_temp)):
    return paired_ncc_file_name, ambig_ncc_file_name

  # Load a list of unpairable hybrid segments

  unpaired_dict = defaultdict(set)
  unpaired_size = None

  for file_path in (unpair_path1, unpair_path2):
    with open_file_r(file_path) as file_obj:
      for line in file_obj:
        contig, start, end = line.split()
        start = int(start)

        if not unpaired_size:
          unpaired_size = int(end) - start

        unpaired_dict[contig].add(start)

  unpaired_offset = unpaired_size/2

  ncc_file_obj = open(paired_ncc_file_name_temp, 'w')
  write_pair = ncc_file_obj.write

  ambig_ncc_file = open(ambig_ncc_file_name_temp, 'w')
  write_ambig = ambig_ncc_file.write

  file_objs = [open_file_r(sam_file) for sam_file in (sam_file1, sam_file2, sam_file3, sam_file4)]
  readlines = [f.readline for f in file_objs]

  n_map = [0,0,0,0]
  n_pairs = 0
  n_ambig = 0
  n_unpaired = 0
  n_unambig = 0
  n_unmapped = 0
  n_hybrid_ambig = 0
  n_unpairable = 0

  # Go through same files and pair based on matching id
  # Write out any multi-position mapings to ambiguous

  # Skip to end of headers

  lines = [_skip_sam_header(readline) for readline in readlines]

  # Process data lines

  ids = [line[:10] for line in lines]

  while '' not in ids:
    _id = max(ids)

    for i in range(4):
      if ids[i] < _id: # Reads may not match as some ends can be removed during clipping
        n_unpaired += 1
        line = readlines[i]()
        lines[i] = line
        ids[i] = line[:10]
        break

    else:
      contacts = [[], [], [], []]
      scores = [[], [], [], []]

      for i in range(4):
        while ids[i] == _id:
          n_map[i] += 1
          line = lines[i]
          data = line.split('\t')
          chromo = data[2]
          strand = '-' if int(data[1]) & 0x10 else '+'
          start = int(data[3])
          seq = data[9]
          end = start + len(seq)

          if strand == '-':
            start, end = end, start  # The sequencing read started from the other end

          if chromo != '*':
            score = int(SCORE_TAG_SEARCH(line).group(1))
            ncc = (chromo, 0, 0, start, end, strand) # Strand info kept because ends can be diffferent for replicate reads, no Re fragment positions, yet
            contacts[i].append((ncc, score))
            scores[i].append(score)

          line = readlines[i]()
          lines[i] = line
          ids[i] = line[:10]

      n_pairs += 1

      if not unique_map: # Resolve some perfect vs non-perfect positional ambiguity
        for a, b in ((0,1), (2,3)):
          if (len(scores[a]) > 1) or (len(scores[b]) > 1):
            n_best_a = scores[a].count(0)# Zero is best/perfect score in end-to-end mode
            n_best_b = scores[b].count(0)

            if n_best_a * n_best_b == 1: # Only one perfect pair
              i = scores[a].index(0)
              j = scores[b].index(0)

              contacts[a] = [contacts[a][i]]
              contacts[b] = [contacts[b][j]]

      # Skip anything identifyably unpairable across genomes
      unpairable = [False, False, False, False]

      for i, end_poss in enumerate(contacts):
        for ncc, score in end_poss:
          chromo = ncc[0]
          start, end = sorted(ncc[3:5])

          j = unpaired_offset * int(start/unpaired_offset)
          n_unpair = 0

          while j < end:
            if j in unpaired_dict[chromo]:
              n_unpair += 1

            j += unpaired_offset

          if n_unpair > max_unpairable:
            unpairable[i] = True
            break

      # If all internal sections of any end unmappable across genomes, then skip

      n0 = len(contacts[0])
      n1 = len(contacts[1])
      n2 = len(contacts[2])
      n3 = len(contacts[3])
      m = max(n0, n1, n2, n3)

      iid = int(_id)

      if (n0+n2) * (n1+n3) == 0: # One end or both ends have no mappings
        n_unmapped += 1
        continue

      elif m == 1: # No multi position mapping

        for end_a, end_b in ((0,2),(2,0),(1,3),(3,1)): # Same reads
          if contacts[end_a] and unpairable[end_a]: # Mapped to one genome but other could be undetactable
            if not contacts[end_b]: # Nothing detected for other genome
              n_unpairable += 1
              break

        else: # Survived the hybrid pairability check

          write_func = write_pair # Send to main file
          pair_scores_a = []
          pair_scores_b = []

          for end_a in (0,2):
            for ncc_a, score_a in contacts[end_a]:
              pair_scores_a.append(score_a)

          for end_b in (1,3):
            for ncc_b, score_b in contacts[end_b]:
              pair_scores_b.append(score_b)

          max_a = max(pair_scores_a)
          max_b = max(pair_scores_b)
          ap = 0
          pairs = []

          for end_a, end_b in ((0,1), (2,3), (0,3), (2,1)): # A:A, B:B, A:B, B:A genome pairings
            for ncc_a, score_a in contacts[end_a]:
              if score_a < max_a:
                continue

              for ncc_b, score_b in contacts[end_b]:
                if score_b < max_b:
                  continue

                ap += 1
                pairs.append((end_a, end_b, ncc_a, ncc_b))

          if ap > 1:
            for end_a, end_b, ncc_a, ncc_b in pairs: # Allow imperfect matches if genome ambigous
              _write_ncc_line(ncc_a, ncc_b, n_pairs, iid, write_func)

            n_hybrid_ambig += 1

          else:
            if (max_a == 0) and (max_b == 0): # Only allow perfect unambigous
              for end_a, end_b, ncc_a, ncc_b in pairs:
                _write_ncc_line(ncc_a, ncc_b, n_pairs, iid, write_func)

              n_unambig += 1

            else:
              n_unpairable += 1

      else:
        n_ambig += 1
        write_func = write_ambig # Send to ambig file

        for end_a, end_b in ((0,1), (2,3), (0,3), (2,1)): # A:A, B:B, A:B, B:A genome pairings
          for ncc_a, score_a in contacts[end_a]:
            for ncc_b, score_b in contacts[end_b]:
              _write_ncc_line(ncc_a, ncc_b, n_pairs, iid, write_func)

  ncc_file_obj.close()

  stats = [('end_1A_alignments', n_map[0]),
           ('end_2A_alignments', n_map[1]),
           ('end_1B_alignments', n_map[2]),
           ('end_2B_alignments', n_map[3]),
           ('unpaired_ends', n_unpaired),
           ('hybrid_unpairable', (n_unpairable, n_pairs)),
           ('unmapped_end', (n_unmapped, n_pairs)),
           ('unique', (n_unambig, n_pairs)),
           ('genome_ambiguous', (n_hybrid_ambig, n_pairs)),
           ('position_ambiguous', (n_ambig, n_pairs)),
           ('total_contacts', n_pairs)]

  log_report('pair', stats)

  move(paired_ncc_file_name_temp, paired_ncc_file_name)
  move(ambig_ncc_file_name_temp, ambig_ncc_file_name)

  return paired_ncc_file_name, ambig_ncc_file_name


def pair_mapped_seqs(sam_file1, sam_file2, file_root, ambig=True, unique_map=False):

  paired_ncc_file_name = tag_file_name(file_root, 'pair', '.ncc')
  ambig_ncc_file_name = tag_file_name(file_root, 'pair_ambig', '.ncc')
  paired_ncc_file_name_temp = paired_ncc_file_name + TEMP_EXT
  ambig_ncc_file_name_temp = ambig_ncc_file_name + TEMP_EXT

  if (INTERRUPTED and os.path.exists(paired_ncc_file_name)
      and os.path.exists(ambig_ncc_file_name)
      and not os.path.exists(paired_ncc_file_name_temp)
      and not os.path.exists(ambig_ncc_file_name_temp)):
    return paired_ncc_file_name, ambig_ncc_file_name

  ncc_file_obj = open(paired_ncc_file_name_temp, 'w')
  write_pair = ncc_file_obj.write

  ambig_ncc_file = open(ambig_ncc_file_name_temp, 'w')
  write_ambig = ambig_ncc_file.write

  file_obj1 = open_file_r(sam_file1)
  file_obj2 = open_file_r(sam_file2)

  readline1 = file_obj1.readline
  readline2 = file_obj2.readline

  n_map1 = 0
  n_map2 = 0
  n_pairs = 0
  n_ambig = 0
  n_unpaired = 0
  n_unambig = 0
  n_unmapped = 0

  # Go through same files and pair based on matching id
  # Write out any to-many mapings to ambiguous

  # Skip to end of headers

  line1 = _skip_sam_header(readline1)
  line2 = _skip_sam_header(readline2)

  # Process data lines

  id1 = line1[:10]
  id2 = line2[:10]

  while id1 and id2:
    if id1 < id2:
      n_unpaired += 1
      line1 = readline1()
      id1 = line1[:10]
      continue

    if id2 < id1:
      n_unpaired += 1
      line2 = readline2()
      id2 = line2[:10]
      continue

    _id = id1

    contact_a = []
    contact_b = []
    scores_a = []
    scores_b = []

    if id1 != id2:
      raise Exception('Pair mismatch')

    while id1 == _id:
      n_map1 += 1
      data_a = line1.split('\t')
      chr_a = data_a[2]
      strand_a = '-' if int(data_a[1]) & 0x10 else '+'
      start_a = int(data_a[3])
      seq_a = data_a[9]
      end_a = start_a + len(seq_a)

      if strand_a == '-':
        start_a, end_a = end_a, start_a  # The sequencing read started from the other end

      if chr_a != '*':
        score_a = int(SCORE_TAG_SEARCH(line1).group(1))
        ncc_a = (chr_a, 0, 0, start_a, end_a, strand_a) # Strand info kept because ends can be diffferent for replicate reads, no Re fragment positions, yet
        contact_a.append((ncc_a, score_a))
        scores_a.append(score_a)

      line1 = readline1()
      id1 = line1[:10]

    while id2 == _id:
      n_map2 += 1
      data_b = line2.split('\t')
      chr_b = data_b[2]
      strand_b = '-' if int(data_b[1]) & 0x10 else '+'
      start_b = int(data_b[3])
      seq_b = data_b[9]
      end_b = start_b + len(seq_b)

      if strand_b == '-':
        start_b, end_b = end_b, start_b  # The sequencing read started from the other end

      if chr_b != '*':
        score_b = int(SCORE_TAG_SEARCH(line2).group(1))
        ncc_b = (chr_b, 0, 0, start_b, end_b, strand_b)
        contact_b.append((ncc_b, score_b))
        scores_b.append(score_b)

      line2 = readline2()
      id2 = line2[:10]

    if not unique_map: # Resolve some perfect vs non-perfect positional ambiguity
      if (len(scores_a) > 1) or (len(scores_b) > 1):
        n_best_a = scores_a.count(0)# Zero is best/perfect score in end-to-end mode
        n_best_b = scores_b.count(0)

        if n_best_a * n_best_b == 1: # Only one perfect pair
          i = scores_a.index(0)
          j = scores_b.index(0)

          contact_a = [contact_a[i]]
          contact_b = [contact_b[j]]

    n_pairs += 1
    n_a = len(contact_a)
    n_b = len(contact_b)
    m = max(n_a, n_b)
    iid = int(_id)

    if n_a * n_b == 0: # Either end not mapped
      n_unmapped += 1
      continue

    elif m == 1: # No multi-position mapping
      n_unambig += 1
      write_func = write_pair # Send to main file

    else:
      n_ambig += 1
      write_func = write_ambig # Send to ambig file

    for ncc_a, score_a in contact_a:
      for ncc_b, score_b in contact_b:
        _write_ncc_line(ncc_a, ncc_b, n_pairs, iid, write_func) # Write all pair combinations where ambigous

  ncc_file_obj.close()

  stats = [('end_1_alignments', n_map1),
           ('end_2_alignments', n_map2),
           ('unpaired_ends', n_unpaired),
           ('unmapped_end', (n_unmapped, n_pairs)),
           ('unique', (n_unambig, n_pairs)),
           ('ambiguous', (n_ambig, n_pairs)),
           ('total_contacts', n_pairs)]

  log_report('pair', stats)

  move(paired_ncc_file_name_temp, paired_ncc_file_name)
  move(ambig_ncc_file_name_temp, ambig_ncc_file_name)

  return paired_ncc_file_name, ambig_ncc_file_name


def map_reads(fastq_file, genome_index, align_exe, num_cpu, ambig, qual_scheme, job):

  sam_file_path = tag_file_name(fastq_file, 'map%d' % job, '.sam')
  sam_file_path_temp = sam_file_path + TEMP_EXT

  if INTERRUPTED and os.path.exists(sam_file_path) and not os.path.exists(sam_file_path_temp):
    return sam_file_path

  patt_1 = re.compile('(\d+) reads; of these:')
  patt_2 = re.compile('(\d+) \(.+\) aligned exactly 1 time')
  patt_3 = re.compile('(\d+) \(.+\) aligned 0 times')
  patt_4 = re.compile('(\d+) \(.+\) aligned >1 times')

  cmd_args = [align_exe,
              '-D', '20', '-R', '3', '-N', '0',  '-L', '20',  '-i', 'S,1,0.50',
              '-x', genome_index,
              '-k', '2',
              '--reorder',
              '-p', str(num_cpu),
              '-U', fastq_file,
              '-S', sam_file_path_temp]

  if qual_scheme == 'phred33':
    cmd_args += ['--phred33']

  elif qual_scheme == 'phred64':
    cmd_args += ['--phred64']

  elif qual_scheme == 'solexa':
    cmd_args += ['--solexa-quals']

  n_reads = 0
  n_uniq = 0
  n_unmap = 0
  n_ambig = 0

  proc = Popen(cmd_args, stdin=PIPE, stderr=PIPE)
  std_out, std_err = proc.communicate()

  if std_err:
    std_err = std_err.decode('ascii')

    if 'Error' in std_err:
      warn(std_err)

    for line in std_err.split('\n'):

      if line.strip():

        match = patt_1.search(line)
        if match:
          n_reads = int(match.group(1))

        match = patt_2.search(line)
        if match:
          n_uniq = int(match.group(1))

        match = patt_3.search(line)
        if match:
          n_unmap = int(match.group(1))

        match = patt_4.search(line)
        if match:
          n_ambig = int(match.group(1))

  stats = [('input_reads',n_reads),
           ('unique',(n_uniq, n_reads)),
           ('ambiguous',(n_ambig, n_reads)),
           ('unmapped',(n_unmap, n_reads))]

  stat_key = 'map_%d' % job
  log_report(stat_key, stats)

  move(sam_file_path_temp, sam_file_path)

  return sam_file_path


def clip_reads(fastq_file, file_root, junct_seq, replaced_seq, qual_scheme, min_qual, is_second=False, min_len=MIN_READ_LEN):
  """
  Clips reads at ligation junctions, removes anbigous ends and discards very short reads
  """

  tag = 'reads2_clipped' if is_second else 'reads1_clipped'
  clipped_file = tag_file_name(file_root, tag, '.fastq')
  clipped_file_temp = clipped_file + TEMP_EXT

  if INTERRUPTED and os.path.exists(clipped_file) and not os.path.exists(clipped_file_temp):
    return clipped_file

  in_file_obj = open_file_r(fastq_file)
  n_rep = len(replaced_seq)
  n_junc = len(junct_seq)

  n_reads = 0
  n_jclip = 0
  n_qclip = 0
  n_short = 0
  mean_len = 0

  zero_ord = QUAL_ZERO_ORDS[qual_scheme]

  out_file_obj = open(clipped_file_temp, 'w', READ_BUFFER)
  write = out_file_obj.write
  readline = in_file_obj.readline

  line1 = readline()
  while line1[0] != '@':
    line1 = readline()

  while line1:
    n_reads += 1
    line2 = readline()[:-1]
    line3 = readline()
    line4 = readline()[:-1]

    if junct_seq in line2:
      n_jclip += 1
      i = line2.index(junct_seq)
      line2 = line2[:i] + replaced_seq
      line4 = line4[:i+n_rep]

    q = 0
    while line2 and line2[-1] == 'N':
      q = 1
      line2 = line2[:-1]
      line4 = line4[:-1]

    while line4 and (ord(line4[-1]) - zero_ord) < min_qual:
      q = 1
      line2 = line2[:-1]
      line4 = line4[:-1]

    n_qclip += q

    if n_junc < n_rep:
      while len(line4) < len(line2):
        line4 += line4[-1]

    n = len(line2)

    if n < min_len:
      n_short += 1

    else:
      mean_len += len(line2)
      line1 = '@%10.10d_%s' % (n_reads, line1[1:])
      write('%s%s\n%s%s\n' % (line1, line2, line3, line4))

    line1 = readline()

  if n_reads:
    mean_len /= float(n_reads)

  stats = [('input_reads',n_reads),
           ('junction_clipped',(n_jclip, n_reads)),
           ('quaility_clipped',(n_qclip, n_reads)),
           ('too_short',(n_short, n_reads)),
           ('mean_length',mean_len)]

  stat_key = 'clip_2' if is_second else 'clip_1'
  log_report(stat_key, stats)

  move(clipped_file_temp, clipped_file)

  return clipped_file


def get_hybrid_unpairable(fasta_files, genome_index_1, genome_index_2, hom_chromo_dict, align_exe, num_cpu, fragment_length=50, fragment_offset=25):

  # Go through each chromosome of each genome and get a FASTA file of all (sensible) fragments
  # Record which fragment positions cannot be sensibly mapped to any position in the other genome
  # must allow for mismatches and poor matches.
  # be aware of soft-masked repeats and N's
  # need only search between homologous chromosomes

  # When used for filtering it is better if both ends are not observable in a genome
  # compare to observing one from each

  # Output is simply a list of unpairable bins

  # Think about a similar analysis where a sequence (not from the genome build)
  # doesn't map to itself but does map to another chromosome

  genome1 = os.path.basename(genome_index_1)
  genome2 = os.path.basename(genome_index_2)

  unpair_file_name = 'unpaired_%s_%s.tsv' % (genome1, genome2)
  unpair_file_path = os.path.join(os.path.dirname(genome_index_1), unpair_file_name)

  if os.path.exists(unpair_file_path):
    msg = 'Found existing hybrid mappability file %s'
    info(msg % unpair_file_path)
    return unpair_file_path

  if not fasta_files:
    msg = 'No nucleotide sequence FASTA files specified for genome %s: required to check unparable sequences in genome %s'
    fatal(msg % (genome_index_1, genome_index_2))

  # Write FASTQ queries, split so that smaller SAM files are made during mapping check
  msg = 'Calculating hybrid mappability of %s sequences in genome %s'
  info(msg % (genome1, genome2))

  fasta_file_names = {}

  k = 0
  n_bases = 0

  for fasta_file in fasta_files:
    with open_file_r(fasta_file) as in_file_obj:
      seq = ''
      contig = None

      for line in in_file_obj:
        if line[0] == '>':
          if seq and contig and (contig in hom_chromo_dict): # Write out prev
            info('  .. fragmenting {:,} bp sequence for contig {}'.format(len(seq), contig))
            fasta_file_name = '_temp_%s_pair_check.fasta' % (contig)
            fasta_file_names[contig] = fasta_file_name
            n = len(seq)
            n_bases += n

            if os.path.exists(fasta_file_name):
              contig = line[1:].split()[0] # Next header
              seq = ''
              continue

            out_file_obj = open(fasta_file_name, 'w')
            write = out_file_obj.write

            i = 0

            while i < n:
              j = i+fragment_length
              sub_seq = seq[i:j]
              a = b = 0

              while sub_seq and sub_seq[0] == 'N':
                a += 1
                sub_seq = sub_seq[1:]

              while sub_seq and sub_seq[-1] == 'N':
                b += 1
                sub_seq = sub_seq[:-1]

              if len(sub_seq) > MIN_READ_LEN:
                fasta_line = '>%s__%d__%d__%d__%d\n%s\n' % (contig, i, j, a, b, sub_seq)
                write(fasta_line)
                k += 1

              i += fragment_offset

            out_file_obj.close()

          contig = line[1:].split()[0] # Next header
          seq = ''

        else:
          seq += line[:-1]

      if seq and contig and (contig in hom_chromo_dict): # Write out last
        info('  .. found sequence for contig %s' % (contig,))
        fasta_file_name = '_temp_%s_pair_check.fasta' % (contig)
        fasta_file_names[contig] = fasta_file_name
        n = len(seq)
        n_bases += n

        if os.path.exists(fasta_file_name):
          contig = line[1:].split()[0] # Next header
          seq = ''
          continue

        out_file_obj = open(fasta_file_name, 'w')
        write = out_file_obj.write

        i = 0

        while i < n:
          j = i+fragment_length
          sub_seq = seq[i:j]
          a = b = 0

          while sub_seq and sub_seq[0] == 'N':
            a += 1
            sub_seq = sub_seq[1:]

          while sub_seq and sub_seq[-1] == 'N':
            b += 1
            sub_seq = sub_seq[:-1]

          if len(sub_seq) > MIN_READ_LEN:
            fasta_line = '>%s__%d__%d__%d__%d\n%s\n' % (contig, i, j, a, b, sub_seq)
            write(fasta_line)
            k += 1

          i += fragment_offset

        out_file_obj.close()

  contigs = sorted(fasta_file_names)
  msg = 'Calculating hybrid mappability for {:,} fragments of {:,} sequences totalling {:,} bp'
  info(msg.format(k, len(contigs), n_bases))

  with open(unpair_file_path, 'w') as unpair_file_obj:
    unpair_write = unpair_file_obj.write

    for i, contig in enumerate(contigs):
      info('  .. contig %s - %d of %d' % (contig, i+1, len(contigs)))

      sam_file_path = tag_file_name(fasta_file_names[contig], '', '.sam')

      cmd_args = [align_exe, '-D', '20', '-R', '3', '-N', '1', # Allow a mismatch
                  '-L', '20', '-i', 'S,1,0.50',
                  '-x', genome_index_2,
                  '-p', str(num_cpu),
                  '-f', '-U', fasta_file_names[contig],
                  '-S', sam_file_path]

      proc = Popen(cmd_args, stdin=PIPE)
      std_out, std_err = proc.communicate()

      sam_file_obj = open(sam_file_path, 'r')

      for line in sam_file_obj:
        if line[0] == '@':
          continue

        sam_data = line.split('\t')
        if not sam_data:
          continue

        sam_bits = int(sam_data[1])

        if sam_bits & 0x4: # Not mapped
          fasta_head = sam_data[0]
          contig, start, end, off1, off2 = fasta_head.split('__')
          unpair_write('%s\t%s\t%s\n' % (contig, start, end))

      sam_file_obj.close()
      os.unlink(sam_file_path)
      os.unlink(fasta_file_names[contig])

  return unpair_file_path


def get_chromo_re_fragments(fasta_file_objs, contig, sequence, re_site, cut_offset, mappability_length=75):

  frag_data = []
  frag_data_append = frag_data.append
  re_site = re_site.replace('^', '')
  re_site = re_site.replace('_', '')

  if 'N' in re_site:
    re_site = re_site.replace('N', '[AGCT]') # Could allow more ambiguity codes, if feeling generous
    frags = re.split(re_site, sequence)
  else:
    frags = sequence.split(re_site)

  n_files = len(fasta_file_objs)

  # Make FASTA file of fragments to assess mappability

  fasta_write = [fo.write for fo in fasta_file_objs]

  step = mappability_length/2
  site_len = len(re_site)
  offset_start = site_len - cut_offset
  frag_start = offset_start
  k = 0 # Indicates which FASTA file to write to

  for frag_seq in frags:
    n = len(frag_seq)

    start = frag_start - offset_start
    end = frag_start + n + cut_offset # Offset is cut site, which can be in subseq removed during split

    if not n:
      frag_start += site_len
      continue

    a = 0
    b = mappability_length

    while b < n:
      sub_seq = frag_seq[a:b]
      m = b-a

      # Deal with "N"s

      i = 0 # Left clip
      j = m # Right clip

      while (i < m) and (sub_seq[i] == 'N'):
        i += 1

      while (j > i) and (sub_seq[j-1] == 'N'):
        j -= 1

      sub_seq = sub_seq[i:j]

      if len(sub_seq) > MIN_READ_LEN:
        fasta_line = '>%s__%s__%s\n%s\n' % (contig, start, a+i, sub_seq)
        fasta_write[k % n_files](fasta_line)
        k += 1

      a += step
      b = a + mappability_length

    if b != n:
      a = n-mappability_length
      b = n

      sub_seq = frag_seq[a:b]
      m = len(sub_seq)

      # Deal with "N"s

      i = 0 # Left clip
      j = m # Right clip

      while (i < m) and (sub_seq[i] == 'N'):
        i += 1

      while (j > i+1) and (sub_seq[j-1] == 'N'):
        j -= 1

      sub_seq = sub_seq[i:j]

      if len(sub_seq) > MIN_READ_LEN:
        fasta_line = '>%s__%s__%s\n%s\n' % (contig, start, a+i, sub_seq)
        fasta_write[k % n_files](fasta_line)
        k += 1

    frag_data_append([start, end-1, 0]) # Mappability set later
    frag_start += n + site_len

  frag_data[-1][1] = len(sequence)-1 # Last position is seq end, not cut site

  return frag_data


def get_re_seq(re_name):

  site = RE_SITES[re_name]
  return site.replace('^', '').replace('_','')


def check_re_frag_file(genome_index, re_name, genome_fastas, align_exe, num_cpu, remap=False):

  dir_name, index_file = os.path.split(genome_index)
  re_site_file = 'RE_frag_%s_%s.txt' % (re_name, index_file)
  re_site_file = os.path.join(dir_name, re_site_file)

  if os.path.exists(re_site_file):
    if remap:
      msg = 'Replacing existing restriction enzyme fragment file %s'
      info(msg % re_site_file)

    else:
      msg = 'Found existing restriction enzyme fragment file %s'
      info(msg % re_site_file)
      return re_site_file

  if not genome_fastas:
    msg = 'Restriction enzyme fragment file %s needs to be built but no genome FASTA files have been specified'
    fatal(msg % re_site_file)

  else:
    msg = 'Creating restriction enzyme fragment file %s'
    info(msg % re_site_file)

  re_site = RE_SITES[re_name]
  seq = get_re_seq(re_name)
  cut_pos = re_site.replace('_', '').index('^')

  # FASTA file of geneome fragments split over multiple files for mappability calculation
  # so that none of the temporary SAM files get too big

  fasta_file_names = []
  fasta_file_objs = []

  for i in range(NUM_MAP_FASTAS):
    fasta_file_name = genome_index + '_temp_%03d_re_frags.fasta' % (i+1,)
    fasta_file_obj = open(fasta_file_name, 'w')
    fasta_file_names.append(fasta_file_name)
    fasta_file_objs.append(fasta_file_obj)

  frag_data = {}
  sources = {}

  msg = 'Calculating %s fragment locations for %s contig %s'
  n_bases = 0

  for fasta_file in genome_fastas:
    file_obj = open_file_r(fasta_file)

    contig = None

    if (len(genome_fastas) < 3) and ('_chr' not in os.path.basename(fasta_file)):
      source = os.path.splitext(os.path.basename(fasta_file))[0].split('.')[0]
    else:
      source = os.path.splitext(fasta_file)[0].split('_')[-1]

    seq = ''

    for fasta_line in file_obj:
      if fasta_line[0] == '>':
        if source and seq: # add previous
          info(msg % (re_name, source, contig))
          n_bases += len(seq)
          frag_data[contig] = get_chromo_re_fragments(fasta_file_objs, contig, seq.upper(), re_site, cut_pos)

        contig = fasta_line[1:].split()[0]
        sources[contig] = source
        seq = ''

      else:
        seq += fasta_line[:-1]

    file_obj.close()

    if source and seq:
      info(msg % (re_name, source, contig))
      n_bases += len(seq)
      frag_data[contig] = get_chromo_re_fragments(fasta_file_objs, contig, seq.upper(), re_site, cut_pos)

  for fasta_file_obj in fasta_file_objs:
    fasta_file_obj.close()

  # Calculate mappability

  contigs = sorted(frag_data)
  frag_read_counts = {}

  for contig in contigs:
    frag_read_counts[contig] = defaultdict(int)

  msg = 'Calculating mappability for {:,} sequences totalling {:,} bp'
  info(msg.format(len(contigs), n_bases))

  for i, fasta_file_name in enumerate(fasta_file_names):
    info('  .. chunk %d of %d' % (i+1, NUM_MAP_FASTAS))

    sam_file_path = tag_file_name(fasta_file_name, '', '.sam')

    cmd_args = [align_exe, '--sensitive',
                '-x', genome_index,
                '-p', str(num_cpu),
                '-f', '-U', fasta_file_name,
                '-S', sam_file_path]  # Include unaligned fragments, naturally

    proc = Popen(cmd_args, stdin=PIPE, stderr=PIPE)
    std_out, std_err = proc.communicate()

    sam_file_obj = open(sam_file_path, 'r')

    for line in sam_file_obj:
      if line[0] == '@':
        continue

      sam_data = line.split('\t')
      if not sam_data:
        continue

      sam_bits = int(sam_data[1])

      if sam_bits & 0x4: # Not mapped
        continue

      fasta_head = sam_data[0]
      contig, start, indent = fasta_head.split('__')

      frag_read_counts[contig][int(start)] += 1

    sam_file_obj.close()
    os.unlink(sam_file_path)
    os.unlink(fasta_file_name)

  # Add mappability

  for contig in contigs:
    for i, (start, end, null) in enumerate(frag_data[contig]):
      frag_data[contig][i][2] = frag_read_counts[contig][start]

  # Write out RE file

  out_file_obj = open(re_site_file, 'w')
  head = 'Source\tContig\tFrag_Start\tFrag_End\tMappability\n'
  out_file_obj.write(head)

  for contig in contigs:
    source = sources[contig]

    for start, end, mappability in frag_data[contig]:
      frag_line = '%s\t%s\t%d\t%d\t%d\n' % (source, contig, start, end, mappability)
      out_file_obj.write(frag_line)

  out_file_obj.close()

  return re_site_file


def read_re_frag_file(file_path):

  frag_dict = {}
  chromo_contigs = defaultdict(set)
  npy_cache = file_path + '.npz'

  if os.path.exists(npy_cache):
    in_data_dict = np.load(npy_cache)

    for key in in_data_dict:
      source, contig = key.split()
      chromo_contigs[source].add(contig)
      frag_dict[contig] = in_data_dict[key]

  else:
    file_obj = open_file_r(file_path)
    file_obj.readline()

    prev_contig = None
    for line in file_obj:
      source, contig, start, end, mappability = line.split()
      source = source.split()[-1]

      if contig != prev_contig:
        chromo_contigs[source].add(contig)
        frag_dict[contig] = []
        append = frag_dict[contig].append

      append((int(start), int(end), int(mappability)))
      prev_contig = contig

    file_obj.close()

    kw_args = {}
    for source in chromo_contigs:
      for contig in chromo_contigs[source]:
        frag_dict[contig] = np.array(frag_dict[contig])
        key = '%s %s' % (source, contig)
        kw_args[key] = frag_dict[contig]

    try:
      np.savez(npy_cache, **kw_args)
    except Exception:
      pass # Not to worry if cache fails: it's only a speedup

  end_dict = {}
  for contig in frag_dict:
    end_dict[contig] = frag_dict[contig][:,1]

  chromo_name_dict = {} # Safely mappable contig names
  for source in chromo_contigs:
    contigs = chromo_contigs[source]

    if len(contigs) == 1:
      contig = contigs.pop()
      chromo_name_dict[contig] = source
      chromo_name_dict[source] = source

      frag_dict[source] = frag_dict[contig]
      end_dict[source] = end_dict[contig]

  return frag_dict, end_dict, chromo_name_dict


def uncompress_gz_file(file_name):

  if file_name.endswith('.gz'):
    in_file_obj = gzip.open(file_name, 'rb')

    file_name = file_name[:-3]
    out_file_obj = open(file_name, 'w')

    for line in in_file_obj:
      out_file_obj.write(line)

    # Faster alternative, given sufficient memory:
    # out_file_obj.write(infile_obj.read())

    in_file_obj.close()
    out_file_obj.close()

  return file_name


def index_genome(base_name, file_names, output_dir, indexer_exe='bowtie2-build',
                 table_size=10, quiet=True, pack=True):

  fasta_files = []
  for file_name in file_names:
    file_name = uncompress_gz_file(file_name)
    fasta_files.append(file_name)

  fasta_file_str = ','.join(fasta_files)

  cmd_args = [indexer_exe, '-f']

  if quiet:
    cmd_args.append('-q')

  if pack:
    cmd_args.append('-p')

  cmd_args += ['-t', str(table_size), fasta_file_str, base_name]

  call(cmd_args, cwd=output_dir)


def get_ligation_junction(re_site):
  """
  AG^CT gives AG + CT = AG + CT
              TC   GA

  ^GATC_ gives ---- + GATC = GATC + GATC
               CTAG   ----

  A^GATC_T gives A---- + GATCT = A + GATC + GATC + T
                 TCTAG   ----A

  A_GATC^T gives AGATC + ----T = A + GATC + GATC + T
                 T----   CTAGA
  """

  if '_' in re_site:
    i, j = sorted([re_site.index('^'), re_site.index('_')])
    seq = re_site[:i] + re_site[i+1:j] + re_site[i+1:j] + re_site[j+1:]

  else: # Blunt
    seq = re_site.replace('^', '')

  return seq


def check_re_site(seq):

  if set(seq) - set('^_GCATN'):
    msg = 'Restriction site specification (%s) may only contain characters "G", "C", "A", "T", "N", "_", or "^"' % seq
    return False, msg

  if '^' not in seq:
    msg = 'Restriction site specification (%s) does not contain a "^" cut point' % seq
    return False, msg

  if seq.count('^') > 1:
    if '_' in seq:
      msg = 'Restriction site specification (%s) may only contain one "^" 5\' cut point' % seq

    else:
      msg = 'Restriction site specification (%s) may only contain one "^" 5\' cut point, use "_" for 3\'' % seq

    return False, msg

  if seq.count('_') > 1:
    msg = 'Restriction site specification (%s) may only contain one "_" 3\' cut point' % seq
    return False, msg

  return True, ''


def check_fastq_file(file_path):

  file_obj = open_file_r(file_path)

  lines = file_obj.readlines(FASTQ_READ_CHUNK)
  lines = [l for l in lines if l.rstrip()]

  for i in range(0, len(lines), 4):
    if not (lines[i][0] == '@') and (lines[i+2][0] == '+'):
      msg = 'File "%s" does not appear to be in FASTQ format'
      return False, msg % file_path

  return True, ''


def is_genome_indexed(file_path, file_exts=('.1.bt2','.1.bt2l')):

  for file_ext in file_exts:
    if os.path.exists(file_path + file_ext):
      return True

  return False


def check_index_file(file_path, sub_files=('.1', '.2', '.3', '.4', '.rev.1', '.rev.2')):

  if os.path.exists(file_path + '.1.bt2l'):
    file_ext = '.bt2l' # All build files should be long
  else:
    file_ext = '.bt2'

  for sub_file in sub_files:
    full_path = file_path + sub_file + file_ext
    is_ok, msg = check_regular_file(full_path)

    if not is_ok:
      return is_ok, 'Genome index error. ' + msg

  return True, ''


def check_file_extension(file_path, file_ext):

  file_ext = file_ext.lower()

  if not file_ext.startswith('.'):
    file_ext = '.' + file_ext

  if not file_path.lower().endswith(file_ext):
    file_path = os.path.splitext(file_path)[0] + file_ext

  return file_path


def check_regular_file(file_path):

  msg = ''

  if not os.path.exists(file_path):
    msg = 'File "%s" does not exist' % file_path
    return False, msg

  if not os.path.isfile(file_path):
    msg = 'Location "%s" is not a regular file' % file_path
    return False, msg

  if os.stat(file_path).st_size == 0:
    msg = 'File "%s" is of zero size ' % file_path
    return False, msg

  if not os.access(file_path, os.R_OK):
    msg = 'File "%s" is not readable' % file_path
    return False, msg

  return True, msg


def _write_log_lines(lines, verbose=VERBOSE):

  if LOG_FILE_PATH:
    with open(LOG_FILE_PATH, 'a') as file_obj:
      file_obj.write('\n'.join(lines) + '\n')

  if verbose:
    for line in lines:
      print(line)


def info(msg, prefix='INFO'):

  line = '%8s : %s' % (prefix, msg)
  _write_log_lines([line])


def warn(msg, prefix='WARNING'):

  line = '%8s : %s' % (prefix, msg)
  _write_log_lines([line])


def fatal(msg, prefix='%s FAILURE' % PROG_NAME):

  lines = ['%8s : %s' % (prefix, msg),
           'Use %s with -h to display command line options' % PROG_NAME]

  _write_log_lines(lines, verbose=True)
  sys.exit(0)


def pair_fastq_files(fastq_paths, pair_tags=('r_1','r_2'), err_msg='Could not pair FASTQ read files.'):

  if len(fastq_paths) != len(set(fastq_paths)):
    msg = '%s Repeat file path present.'
    fatal(msg % (err_msg))

  t1, t2 = pair_tags

  paths_1 = []
  paths_2 = []

  for path in fastq_paths:
    dir_name, base_name = os.path.split(path)

    if (t1 in base_name) and (t2 in base_name):
      msg = '%s Tags "%s" and "%s" are ambiguous in file %s'
      fatal(msg % (err_msg, t1, t2, base_name))

    elif t1 in base_name:
      paths_1.append((path, dir_name, base_name))

    elif t2 in base_name:
      paths_2.append((path, dir_name, base_name))

    else:
      msg = '%s File name %s does not contain tag "%s" or "%s"'
      fatal(msg % (err_msg, base_name, t1, t2))

  n1 = len(paths_1)
  n2 = len(paths_2)

  if n1 != n2:
    msg = '%s Number of read 1 (%d) and read 2 (%d) files do not match'
    fatal(msg % (err_msg, n1, n2))

  fastq_paths_r1 = []
  fastq_paths_r2 = []

  for path_1, dir_name_1, file_1 in paths_1:
    seek_file = file_1.replace(t1, t2)
    found = []

    for path_2, dir_name_2, file_2 in paths_2:
      if dir_name_1 != dir_name_2:
        continue

      if file_2 == seek_file:
        found.append(path_2)

    if len(found) == 0:
      # No need to check unpaired read 2 files as these always result in an isolated read 1
      msg = '%s No read 2 file "%s" found to pair with %s'
      fatal(msg % (err_msg, seek_file, path_1))

    else:
      # Repeat Read 2 files not possible as repeats checked earlier
      fastq_paths_r1.append(path_1)
      fastq_paths_r2.append(found[0])

  return fastq_paths_r1, fastq_paths_r2


def is_interrupted_job():

  if STAT_FILE_PATH and os.path.exists(STAT_FILE_PATH):
    with open(STAT_FILE_PATH) as file_obj:
      stat_dict = json.load(file_obj)

    if 'final' in stat_dict:
      del stat_dict['final']

      with open(STAT_FILE_PATH, 'w') as file_obj: # Overwrite
        json.dump(stat_dict, file_obj)

      if os.path.exists(LOG_FILE_PATH):
        os.remove(LOG_FILE_PATH)

      return False

    else:
      info('* * Resume interrupted job * * ')
      return True

  else:
    if os.path.exists(LOG_FILE_PATH):
      os.remove(LOG_FILE_PATH)

    return False


def get_fastq_qual_scheme(file_path):
  """
  Guess the quality scoring scheme for a FASTQ file
  """

  file_obj = open_file_r(file_path)

  lines = file_obj.readlines(FASTQ_READ_CHUNK)

  while lines[0][0] != '@': # Just in case of headers or other nonsense
    lines.pop(0)

  n_reads = 0
  quals = set()

  n = len(lines)
  n_reads += n // 4

  for i in range(n_reads):
    line_quals = lines[i*4+3][:-1]
    quals.update(set(line_quals))

  file_obj.close()

  quals = [ord(x) for x in quals]
  min_qual = min(quals)
  max_qual = max(quals)

  if min_qual < 33:
    scheme = 'integer'

  elif (max_qual < 75) and (min_qual < 59): # Sanger, Illumina 1.8+
    scheme = 'phred33'

  elif (max_qual > 74):
    if min_qual < 64:
      scheme = 'solexa'
    else:
      scheme = 'phred64'

  else:
    warn('FASTQ quality scheme could not be determined. Assuming Phred+33 (Illumina 1.8+)')
    scheme = 'phred33'

  return scheme


def write_report(report_file, second_genome=None):

  with open(STAT_FILE_PATH) as file_obj:
    stat_dict = json.load(file_obj)

  general_stats = stat_dict['general']
  clip_stats1 = stat_dict['clip_1']
  clip_stats2 = stat_dict['clip_2']
  sam_stats1 = stat_dict['map_1']
  sam_stats2 = stat_dict['map_2']

  if second_genome:
    sam_stats3 = stat_dict['map_3']
    sam_stats4 = stat_dict['map_4']

  pair_stats = stat_dict['pair']
  filter_stats = stat_dict['filter']
  redundancy_stats = stat_dict.get('dup')
  promiscuity_stats = stat_dict.get('promsic')
  final_stats = stat_dict['final']
  frag_sizes = stat_dict['frag_sizes']

  def format_list(d):

    l = []
    for k, v in d:
      c1 = k.replace('_', ' ')
      c1 = c1[0].upper() + c1[1:]

      if isinstance(v, (tuple, list)):
        v, n = v
        percent = 100.0 * v/float(n or 1)
        c2 = '{:,}'.format(v)
        c3 = '{:.2f}%'.format(percent)

      elif isinstance(v, int):
        c2 = '{:,}'.format(v)
        c3 = None

      elif isinstance(v, float):
        c2 = '{:.3f}'.format(v)
        c3 = None

      else:
        c2 = v
        c3 = None

      l.append([c1, c2, c3])

    return l

  def pie_values(data, names):

    sizes = []
    names_2 = []

    data = dict(data)

    for key in names:
      val = data[key]

      if isinstance(val, (tuple, list)):
        val, n = val
      sizes.append(val)
      name = key.replace('_', ' ')
      name = name[0].upper() + name[1:]
      names_2.append(name)

    if sizes and sum(sizes):
      return sizes, names_2

    else:
      return [100], ['No data']

  svg_doc = SvgDocument()

  head_color = 'black'
  main_color = '#006464'
  head_size = 20
  head_height = 32
  head_pad = 7
  main_pad = 16
  row_size = 16
  table_width = 800
  chart_height = 7 * row_size
  text_anchors = ['start','end','end']

  x = main_pad
  x1 = x + table_width/2 + main_pad

  y = 0
  y += head_height

  title = '%s version %s report' % (PROG_NAME, VERSION)
  svg_doc.text(title, (x,y), head_size, bold=True, color=main_color)
  y += head_height

  y += head_pad
  svg_doc.text('Input parameters', (x,y), head_size, bold=True, color=head_color)
  y += head_pad

  data = format_list(general_stats)
  tw, th = svg_doc.table(x, y, table_width, data, False, col_formats=None, size=row_size, main_color=main_color)
  y += th

  y += head_height
  svg_doc.text('Clipping reads 1', (x,y), head_size, bold=True, color=head_color)

  y += head_pad
  data = format_list(clip_stats1)
  tw, th = svg_doc.table(x, y, table_width/2, data, False, text_anchors, col_formats=None, size=row_size, main_color=main_color)
  y += th

  y += head_height
  svg_doc.text('Clipping reads 2', (x,y), head_size, bold=True, color=head_color)

  y += head_pad
  data = format_list(clip_stats2)
  tw, th = svg_doc.table(x, y, table_width/2, data, False, text_anchors, col_formats=None, size=row_size, main_color=main_color)
  y += th

  y += head_height
  svg_doc.text('Genome alignment reads 1', (x,y), head_size, bold=True, color=head_color)

  y += head_pad
  data = format_list(sam_stats1)
  tw, th = svg_doc.table(x, y, table_width/2, data, False, text_anchors, col_formats=None, size=row_size, main_color=main_color)
  names = ['unique','ambiguous','unmapped']
  colors = ['#80C0FF','#FFFF00','#FF0000']
  vals, names = pie_values(sam_stats1, names)
  svg_doc.pie_chart(x1, y, chart_height, vals, names, colors, line_color='#808080')
  y += th

  y += head_height
  svg_doc.text('Genome alignment reads 2', (x,y), head_size, bold=True, color=head_color)

  y += head_pad
  data = format_list(sam_stats2)
  tw, th = svg_doc.table(x, y, table_width/2, data, False, text_anchors, col_formats=None, size=row_size, main_color=main_color)
  names = ['unique','ambiguous','unmapped']
  colors = ['#80C0FF','#FFFF00','#FF0000']
  vals, names = pie_values(sam_stats2, names)
  svg_doc.pie_chart(x1, y, chart_height, vals, names, colors, line_color='#808080')
  y += th

  if second_genome:
    y += head_height
    svg_doc.text('Genome alignment 2 reads 1', (x,y), head_size, bold=True, color=head_color)

    y += head_pad
    data = format_list(sam_stats3)
    tw, th = svg_doc.table(x, y, table_width/2, data, False, text_anchors, col_formats=None, size=row_size, main_color=main_color)
    names = ['unique','ambiguous','unmapped']
    colors = ['#80C0FF','#FFFF00','#FF0000']
    vals, names = pie_values(sam_stats1, names)
    svg_doc.pie_chart(x1, y, chart_height, vals, names, colors, line_color='#808080')
    y += th

    y += head_height
    svg_doc.text('Genome alignment 2 reads 2', (x,y), head_size, bold=True, color=head_color)

    y += head_pad
    data = format_list(sam_stats4)
    tw, th = svg_doc.table(x, y, table_width/2, data, False, text_anchors, col_formats=None, size=row_size, main_color=main_color)
    names = ['unique','ambiguous','unmapped']
    colors = ['#80C0FF','#FFFF00','#FF0000']
    vals, names = pie_values(sam_stats2, names)
    svg_doc.pie_chart(x1, y, chart_height, vals, names, colors, line_color='#808080')
    y += th

  y += head_height
  svg_doc.text('Pairing reads', (x,y), head_size, bold=True, color=head_color)

  y += head_pad
  data = format_list(pair_stats)
  tw, th = svg_doc.table(x, y, table_width/2, data, False, text_anchors, col_formats=None, size=row_size, main_color=main_color)

  if second_genome:
    names = ['unique','genome_ambiguous','position_ambiguous','unmapped_end']
    colors = ['#80C0FF','#C0C0C0','#FFFF00','#FF0000']

  else:
    names = ['unique','ambiguous','unmapped_end']
    colors = ['#80C0FF','#FFFF00','#FF0000']

  vals, names = pie_values(pair_stats, names)
  svg_doc.pie_chart(x1, y, chart_height, vals, names, colors, line_color='#808080')

  y += th

  y += head_height
  svg_doc.text('Filtering pairs', (x,y), head_size, bold=True, color=head_color)

  y += head_pad
  data = format_list(filter_stats)
  tw, th = svg_doc.table(x, y, table_width/2, data, False, text_anchors, col_formats=None, size=row_size, main_color=main_color)
  names = ['accepted','internal_re1','adjacent_re1','circular_re1', 'overhang_re1',
           'too_close','too_small','too_big',
           'internal_re2','no_end_re2', 'unknown_contig']
  colors = ['#80C0FF','#FFFF00','#FF0000','#C0C0C0','#FF0080',
            '#00FFFF','#00C0C0','#008080',
            '#40FF40','#008000','#404040']
  vals, names = pie_values(filter_stats, names)
  svg_doc.pie_chart(x1, y, chart_height, vals, names, colors, line_color='#808080', small_val=0.01)

  hist, edges = frag_sizes
  data_set = [(int(edges[i]), int(val)) for i, val in enumerate(hist)]
  svg_doc.graph(x1, y+chart_height, table_width/2, th-chart_height-40, [data_set], 'Size', 'Count',
                names=None, colors=None,  graph_type='line',
                symbols=None, line_widths=None, symbol_sizes=None,
                legend=False, title=None, plot_offset=(80, 20))
  y += th

  if redundancy_stats is not None:
    y += head_height
    svg_doc.text('Duplicate removal', (x,y), head_size, bold=True, color=head_color)

    y += head_pad
    data = format_list(redundancy_stats)
    tw, th = svg_doc.table(x, y, table_width/2, data, False, text_anchors, col_formats=None, size=row_size, main_color=main_color)
    names = ['redundant','unique',]
    colors = ['#80C0FF','#FF0000']
    vals, names = pie_values(redundancy_stats, names)
    svg_doc.pie_chart(x1, y, chart_height, vals, names, colors, line_color='#808080')
    y += th

  if promiscuity_stats is not None:
    y += head_height
    svg_doc.text('Promiscuous pair removal', (x,y), head_size, bold=True, color=head_color)

    y += head_pad
    data = format_list(promiscuity_stats)
    tw, th = svg_doc.table(x, y, table_width/2, data, False, text_anchors, col_formats=None, size=row_size, main_color=main_color)
    names = ['clean','resolved','promiscuous']
    colors = ['#80C0FF','#FFFF00','#FF0000']
    vals, names = pie_values(promiscuity_stats, names)
    svg_doc.pie_chart(x1, y, chart_height, vals, names, colors, line_color='#808080')
    y += th

  y += head_height
  svg_doc.text('Final output', (x,y), head_size, bold=True, color=head_color)

  y += head_pad
  data = format_list(final_stats)
  tw, th = svg_doc.table(x, y, table_width/2, data, False, text_anchors, col_formats=None, size=row_size, main_color=main_color)
  names = ['trans','cis_far','cis_near','cis_homolog']
  colors = ['#80C0FF','#FFFF00','#FF0000','#C0C0C0']
  vals, names = pie_values(final_stats, names)
  svg_doc.pie_chart(x1, y, chart_height, vals, names, colors, line_color='#808080')
  y += th

  height = y + main_pad
  width = table_width + 3 * main_pad

  svg_doc.write_file(report_file, width, height)


def log_report(section, data_pairs, extra_dict=None):

  def format_val(val):
    if isinstance(val, int):
      v = '{:,}'.format(val)

    elif isinstance(val, float):
      v = '{:.3f}'.format(val)

    else:
      v = '{}'.format(val)

    return v

  title = STAT_SECTIONS[section]
  lines = [title]

  for key, val in data_pairs:

    if isinstance(val, (tuple, list)):
      val, n = val
      percent = 100.0 * val/float(n or 1.0)
      val = format_val(val)
      info_line = '  {} : {} ({:.2f}%)'.format(key, val, percent)

    else:
      val = format_val(val)
      info_line = '  {} : {}'.format(key, val)

    lines.append(info_line)

  for line in lines:
    info(line)

  if os.path.exists(STAT_FILE_PATH):
    with open(STAT_FILE_PATH) as file_obj:
      stat_dict = json.load(file_obj)

  else:
    stat_dict = {}

  stat_dict[section] = list(data_pairs)

  if extra_dict is not None:
    stat_dict.update(extra_dict)

  with open(STAT_FILE_PATH, 'w') as file_obj: # Overwrite
    json.dump(stat_dict, file_obj)


def read_homologous_chromos(file_path):

  hom_chromo_dict = {}
  haplotype_dict = {}
  chromo_name_dict = {}
  err_msg1 = 'Homologous chromosome file must have 3 whitespace-separated columns for: first parental chromo/contig ID, second parental chromo/contig ID, unified name (e.g. "chrX")'
  err_msg2 = 'Homologous chromosome file contains repeat chromsome/contig name "%s"'

  with open_file_r(file_path) as file_obj:

    for line in file_obj:
      line = line.strip()

      if line[0] == '#':
        continue

      data = line.split()
      if len(data) != 3:
        raise Exception(err_msg1)

      chr_a, chr_b, name = data

      if chr_a in hom_chromo_dict:
        raise Exception(err_msg2 % chr_a)

      if chr_b in hom_chromo_dict:
        raise Exception(err_msg2 % chr_b)

      name_a = name + '.a'
      name_b = name + '.b'

      haplotype_dict[chr_a] = 1
      haplotype_dict[chr_b] = 2

      haplotype_dict[name_a] = 1
      haplotype_dict[name_b] = 2

      hom_chromo_dict[chr_a] = chr_b
      hom_chromo_dict[chr_b] = chr_a

      hom_chromo_dict[name_a] = name_b
      hom_chromo_dict[name_b] = name_a

      chromo_name_dict[chr_a] = name_a
      chromo_name_dict[chr_b] = name_b

  return haplotype_dict, hom_chromo_dict, chromo_name_dict


def nuc_process(fastq_paths, genome_index, genome_index2, re1, re2=None, sizes=(300,800), min_rep=2, num_cpu=1, num_copies=1,
                ambig=True, unique_map=False, homo_chromo=None, out_file=None, ambig_file=None, report_file=None,
                align_exe=None, qual_scheme=None, min_qual=30, g_fastas=None, g_fastas2=None, is_pop_data=False, remap=False, reindex=False,
                keep_files=True, lig_junc=None, zip_files=True, sam_format=True, verbose=True):
  """
  Main function for command-line operation
  """

  genome_indices = [genome_index]
  if genome_index2:
    if not homo_chromo:
      fatal('A homologous chromosome file (-hc option) must be specified for dual genome indices')

    genome_indices.append(genome_index2)

  # Files can be empty if just re-indexing etc...
  if not fastq_paths and not (reindex and g_fastas):
    fatal('No FASTQ files specified')

  # Check FASTQ files : two present, exist, right format

  if fastq_paths:
    if len(fastq_paths) == 1:
      fatal('Only one FASTQ file specified (exactly two required)')

    if len(fastq_paths) > 2:
      fatal('More than two FASTQ files specified (exactly two required)')

  # Ambiguous data is only for single-cell use
  if is_pop_data and ambig:
    fatal('ambiguous option is incompatible with population Hi-C data')

  if not (0 <= min_qual <= 40):
    fatal('Miniumum FASTQ quality score must be in the range 0-40 (%d specified).' % min_qual)

  for file_path in fastq_paths:
    is_ok, msg = check_regular_file(file_path)

    if not is_ok:
      fatal(msg)

    is_ok, msg = check_fastq_file(file_path)

    if not is_ok:
      fatal(msg)

  # Check any homologous chromosome file
  if homo_chromo:
    is_ok, msg = check_regular_file(homo_chromo)

    if not is_ok:
      fatal(msg)

    _, hom_chromo_dict, chromo_name_dict = read_homologous_chromos(homo_chromo)

  else:
    chromo_name_dict = {}
    hom_chromo_dict = {}

  # Check genome index file, if not being rebuilt
  if is_genome_indexed(genome_index):
    is_ok, msg = check_index_file(genome_index)

    if not is_ok:
      fatal(msg)

  if genome_index2:
    if is_genome_indexed(genome_index2):
      is_ok, msg = check_index_file(genome_index2)

      if not is_ok:
        fatal(msg)

  # Check aligner
  if not align_exe:
    try: # Python version >= 3.3
      align_exe = shutil.which('bowtie2')

    except AttributeError:
      cmd_args = ['which', 'bowtie2']
      proc = Popen(cmd_args, stdin=PIPE, stdout=PIPE)
      align_exe, std_err_data = proc.communicate()
      align_exe = align_exe.strip()

  if not align_exe:
    msg = 'Aligner bowtie2 could not be found'
    fatal(msg)

  elif not os.path.exists(align_exe):
    msg = 'Aligner executable path "%s" not found' % align_exe
    fatal(msg)

  else:
    is_ok, msg = check_regular_file(align_exe)

    if not is_ok:
      fatal(msg)

  # Check reindexing input present
  if reindex and not g_fastas:
    fatal('No genome FASTA files specified for re-indexing')

  # Must have a preconstructed genome index or FASTA files to build one from
  if not (genome_index or g_fastas):
    fatal('No genome index file of genome FASTA sequences specified')

  # Check genome FASTA files
  if g_fastas:
    for fata_path in g_fastas:
      is_ok, msg = check_regular_file(fata_path)

      if not is_ok:
        fatal(msg)

  if g_fastas2:
    for fata_path in g_fastas:
      is_ok, msg = check_regular_file(fata_path)

      if not is_ok:
        fatal(msg)

  # Check restriction enzymes
  if re1 not in RE_SITES:
    msg = 'Restriction enzyme "%s" not known. Available: %s.' % (re1, ', '.join(sorted(RE_SITES)))
    msg += ' %s can be edited to add further cut-site defintions.' % RE_CONF_FILE
    fatal(msg)

  is_ok, msg = check_re_site(RE_SITES[re1])
  if not is_ok:
    fatal(msg)

  re1Seq = get_re_seq(re1)

  # FIXME: re2Seq is never used
  if re2:
    re2Seq = get_re_seq(re2)  # noqa: F841

    if re2 not in RE_SITES:
      msg = 'Restriction enzyme "%s" not known. Available: %s.' % (re2, ', '.join(sorted(RE_SITES)))
      msg += ' Enzymes.conf can be edited to add further cut-site defintions.'
      fatal(msg)

    is_ok, msg = check_re_site(RE_SITES[re2])
    if not is_ok:
      fatal(msg)

  else:
    re2Seq = None  # noqa: F841

  # Check FASTQ quality scheme
  if qual_scheme:
    if qual_scheme not in QUAL_SCHEMES:
      msg = 'FASTQ quality scheme "%s" not known. Available: %s.' % (qual_scheme, ', '.join(sorted(QUAL_SCHEMES)))
      msg += ' Scheme will be deduced automatically if not specified.'
      fatal(msg)

  elif fastq_paths:
    qual_scheme = get_fastq_qual_scheme(fastq_paths[0])

  # Check ligation junction sequence
  if lig_junc:
    lig_junc = lig_junc.upper()

    if set(lig_junc) - set('AGCTN'):
      msg = 'Ligation junction sequence may only contain letters G, C, A, T and N'
      fatal(msg)

  else:
    lig_junc = get_ligation_junction(re1Seq)

  # Get base file name for output
  for file_path in (out_file, ambig_file, report_file):
    if file_path:
      file_root = os.path.splitext(file_path)[0]
      break

  else:
    file_paths = list(map(partial(strip_ext, ext=".gz"), fastq_paths))

    merged_path = merge_file_names(file_paths[0], file_paths[1])
    file_root = os.path.splitext(merged_path)[0]

  intermed_dir = file_root + '_nuc_processing_files'
  intermed_file_root = os.path.join(intermed_dir, os.path.basename(file_root))
  if not os.path.exists(intermed_dir):
    os.mkdir(intermed_dir)

  # Check and set output files
  if out_file:
    out_file = check_file_extension(out_file, '.ncc')
  else:
    out_file = file_root + '.ncc'

  if ambig:
    if ambig_file:
      ambig_file = check_file_extension(ambig_file, '.ncc')
    else:
      ambig_file = file_root + '_ambig.ncc'

  else:
    ambig_file = None

  if report_file:
    report_file = check_file_extension(report_file, '.svg')
  else:
    report_file = file_root + '_report.svg'

  global LOG_FILE_PATH
  global STAT_FILE_PATH
  global VERBOSE
  global INTERRUPTED

  LOG_FILE_PATH = file_root + '.log'
  STAT_FILE_PATH = file_root + '_stats.json'
  VERBOSE = bool(verbose)
  INTERRUPTED = is_interrupted_job()

  # Check for no upper fragment size limit
  if len(sizes) == 1:
    size_str = '%d - unlimited' % tuple(sizes)
    sizes = (sizes[0], int(1e16))
  else:
    size_str = '%d - %d' % tuple(sizes)

  general_stats = [
    ('Input FASTQ file 1', fastq_paths[0]),
    ('Input FASTQ file 2', fastq_paths[1]),
    ('Homolog chrom. file', homo_chromo or 'None'),
    ('Output file (main)',os.path.abspath(out_file)),
    ('Output file (ambiguous)',os.path.abspath(ambig_file) if ambig_file else 'None'),
    ('Report file',os.path.abspath(report_file)),
    ('Intermediate file directory',os.path.abspath(intermed_dir)),
    ('Aligner executable',align_exe),
    ('Genome indices',', '.join(genome_indices)),
    ('FASTQ quality scheme',qual_scheme),
    ('Minimum 3\' FASTQ quaility',min_qual),
    ('RE1 site',re1),
    ('RE2 site',(re2 or 'None')),
    ('Ligation junction',lig_junc),
    ('Molecule size range',size_str),
    ('Min sequencing repeats',min_rep),
    ('Parallel CPU cores',num_cpu),
    ('Keep intermediate files?', 'Yes' if keep_files else 'No'),
    ('Input is single-cell Hi-C?', 'No' if is_pop_data else 'Yes'),
    ('Strict mapping only?', 'Yes' if unique_map else 'No'),
    ('SAM output?', 'Yes' if sam_format else 'No'),
    ('GZIP output?', 'Yes' if zip_files else 'No'),
  ]

  log_report('general', general_stats)

  # Check if genome needs indexing

  indexer_exe = os.path.join(os.path.dirname(align_exe), 'bowtie2-build')

  if g_fastas and (reindex or not is_genome_indexed(genome_index)):
    warn('Indexing genome, this may take some time...')
    output_dir, base_name = os.path.split(genome_index)
    index_genome(base_name, g_fastas, output_dir or '.', indexer_exe) # Latest version of bowtie2 can do parallel index builds (--threads)

  if g_fastas2 and genome_index2 and (reindex or not is_genome_indexed(genome_index2)):
    warn('Indexing secondary genome, this may take some time...')
    output_dir, base_name = os.path.split(genome_index2)
    index_genome(base_name, g_fastas2, output_dir or '.', indexer_exe)

  # Create RE fragments file if not present

  re1_files = [check_re_frag_file(genome_index, re1, g_fastas, align_exe, num_cpu, remap=remap)]
  re2_files = []

  if re2:
    re2_files.append(check_re_frag_file(genome_index, re2, g_fastas, align_exe, num_cpu, remap=remap))

  if genome_index2:
    re1_files.append(check_re_frag_file(genome_index2, re1, g_fastas2, align_exe, num_cpu, remap=remap))

    if re2:
      re2_files.append(check_re_frag_file(genome_index2, re2, g_fastas2, align_exe, num_cpu, remap=remap))

    unpair_path1 = get_hybrid_unpairable(g_fastas, genome_index, genome_index2, hom_chromo_dict, align_exe, num_cpu)
    unpair_path2 = get_hybrid_unpairable(g_fastas2, genome_index2, genome_index, hom_chromo_dict, align_exe, num_cpu)

  # Clip read seqs at any sequenced ligation junctions
  info('Clipping FASTQ reads...')
  clipped_file1 = clip_reads(fastq_paths[0], intermed_file_root, lig_junc, re1Seq, qual_scheme, min_qual)
  clipped_file2 = clip_reads(fastq_paths[1], intermed_file_root, lig_junc, re1Seq, qual_scheme, min_qual, is_second=True)

  # Do the main genome mapping
  info('Mapping FASTQ reads...')
  sam_file1 = map_reads(clipped_file1, genome_index, align_exe, num_cpu, ambig, qual_scheme, 1)
  sam_file2 = map_reads(clipped_file2, genome_index, align_exe, num_cpu, ambig, qual_scheme, 2)

  if genome_index2:
    info('Mapping FASTQ reads to second genome...')
    sam_file3 = map_reads(clipped_file1, genome_index2, align_exe, num_cpu, ambig, qual_scheme, 3)
    sam_file4 = map_reads(clipped_file2, genome_index2, align_exe, num_cpu, ambig, qual_scheme, 4)

  if not keep_files:
    os.unlink(clipped_file1)
    os.unlink(clipped_file2)

  info('Pairing FASTQ reads...')

  if genome_index2:
    paired_ncc_file, ambig_paired_ncc_file = pair_mapped_hybrid_seqs(sam_file1, sam_file2, sam_file3, sam_file4, unpair_path1, unpair_path2, intermed_file_root, ambig, unique_map)
  else:
    paired_ncc_file, ambig_paired_ncc_file = pair_mapped_seqs(sam_file1, sam_file2, intermed_file_root, ambig, unique_map)

  # Write SAM, if requested
  if sam_format:
    write_sam_file(paired_ncc_file, sam_file1, sam_file2)

    if ambig:
      write_sam_file(ambig_paired_ncc_file, sam_file1, sam_file2)

  # Filtering
  info('Filtering mapped sequences...')

  filter_output = filter_pairs(paired_ncc_file, re1_files, re2_files, sizes, keep_files, zip_files, chromo_name_dict)
  filter_ncc_file, fail_file_names = filter_output

  # Write SAM, if requested
  if sam_format:
    write_sam_file(filter_ncc_file, sam_file1, sam_file2)

    for key in fail_file_names:
      write_sam_file(fail_file_names[key], sam_file1, sam_file2)

  if ambig:
    info('Filtering ambiguously mapped sequences...')
    filter_output = filter_pairs(ambig_paired_ncc_file, re1_files, re2_files, sizes, keep_files, zip_files, chromo_name_dict, ambig=True)
    ambig_filter_ncc_file, ambig_fail_file_names = filter_output

    # Write SAM, if requested
    if sam_format:
      write_sam_file(ambig_filter_ncc_file, sam_file1, sam_file2)

      for key in fail_file_names:
        write_sam_file(ambig_fail_file_names[key], sam_file1, sam_file2)

  # Redundancy reduction and promiscuity checks not relevant for population data
  # where multiple observations from and between the same RE fragments are expected

  if is_pop_data:
    shutil.copyfile(filter_ncc_file, out_file)

    if sam_format:
      write_sam_file(filter_ncc_file, sam_file1, sam_file2)

    if keep_files:
      if zip_files:
        filter_ncc_file = compress_file(filter_ncc_file)

    else:
      os.unlink(filter_ncc_file)

  else:
    # Merge duplicates
    info('Removing duplicate contacts...')
    nr_ncc_file = remove_redundancy(filter_ncc_file, min_rep, keep_files, zip_files)

    if sam_format:
      write_sam_file(nr_ncc_file, sam_file1, sam_file2)

    if ambig:
      info('Filtering duplicate ambiguous contacts...')
      ambig_nr_ncc_file = remove_redundancy(ambig_filter_ncc_file, min_rep, keep_files, zip_files, ambig=True)

      if sam_format:
        write_sam_file(ambig_nr_ncc_file, sam_file1, sam_file2)

    # Remove promiscuous ends
    info('Removing pairs with promiscuous ends...')
    clean_ncc_file = remove_promiscuous(nr_ncc_file, num_copies, keep_files, zip_files)

    if sam_format:
      write_sam_file(clean_ncc_file, sam_file1, sam_file2)

    shutil.copyfile(clean_ncc_file, out_file)

    if ambig:
      info('Removing ambiguous pairs with promiscuous ends...')
      ambig_clean_ncc_file = remove_promiscuous(ambig_nr_ncc_file, num_copies, keep_files, zip_files, ambig=True)

      if sam_format:
        write_sam_file(ambig_clean_ncc_file, sam_file1, sam_file2)

      shutil.copyfile(ambig_clean_ncc_file, ambig_file)

    if keep_files:
      if zip_files:
        clean_ncc_file = compress_file(clean_ncc_file)
        if ambig:
          ambig_clean_ncc_file = compress_file(ambig_clean_ncc_file)

    else:
      os.unlink(clean_ncc_file)
      if ambig:
        os.unlink(ambig_clean_ncc_file)

  if not keep_files:
    os.unlink(sam_file1)
    os.unlink(sam_file2)

    if genome_index2:
      os.unlink(sam_file3)
      os.unlink(sam_file4)

  final_stats = get_ncc_stats(out_file, hom_chromo_dict)
  log_report('final', final_stats)

  write_report(report_file, genome_index2)

  n_contacts = final_stats[0][1]

  if n_contacts > 1:
    nuc_contact_map(out_file, '_contact_map')

  info('Nuc Process all done.')

  return out_file


def main(argv=None):
  from argparse import ArgumentParser

  if argv is None:
    argv = sys.argv[1:]

  epilog = 'Note %s can be edited to add further restriction enzyme cut-site defintions. ' % RE_CONF_FILE
  epilog += 'For further help email tjs23@cam.ac.uk or wb104@cam.ac.uk'

  res = sorted(RE_SITES)
  avail_re = 'Available: ' + ', '.join(res)
  avail_quals = 'Available: ' + ', '.join(QUAL_SCHEMES)

  arg_parse = ArgumentParser(prog=PROG_NAME, description=DESCRIPTION,
                             epilog=epilog, prefix_chars='-', add_help=True)

  arg_parse.add_argument('-i', nargs='+', metavar='FASTQ_FILE',
                         help='Input paired-read FASTQ files to process. Accepts wildcards that match paired files. If more than two files are input, processing will be run in batch mode using the same parameters.')

  arg_parse.add_argument('-g', metavar='GENOME_FILE',
                         help='Location of genome index files to map sequence reads to without any file extensions like ".1.b2" etc. A new index will be created with the name if the index is missing and genome FASTA files are specified')

  arg_parse.add_argument('-g2', default=None, metavar='GENOME_FILE_2',
                         help='Location of secondary genome index files for hybrid genomes. A new index will be created with the name if the index is missing and genome FASTA files are specified')

  arg_parse.add_argument('-re1', default='MboI', choices=res, metavar='ENZYME',
                         help='Primary restriction enzyme (for ligation junctions). Default: MboI. ' + avail_re)

  arg_parse.add_argument('-re2', choices=res, metavar='ENZYME',
                         help='Secondary restriction enzyme (if used). ' + avail_re)

  arg_parse.add_argument('-s', default='150-1500', metavar='SIZE_RANGE',
                         help='Allowed range of sequenced molecule sizes, e.g. "150-1000", "100,800" or "200" (no maximum)')

  arg_parse.add_argument('-n', default=1, metavar='CPU_COUNT',
                         type=int, help='Number of CPU cores to use in parallel')

  arg_parse.add_argument('-r', default=2, metavar='COUNT',
                         type=int, help='Minimum number of sequencing repeats required to support a contact')

  arg_parse.add_argument('-o', metavar='NCC_FILE',
                         help='Optional output name for NCC format chromosome contact file. This option will be ignored if more than two paired FASTA files are input (i.e. for batch mode); automated naming will be used instead.')

  arg_parse.add_argument('-oa', metavar='NCC_FILE',
                         help='Optional output name for ambiguous contact NCC file. This option will be ignored if more than two paired FASTA files are input (i.e. for batch mode); automated naming will be used instead.')

  arg_parse.add_argument('-or', metavar='REPORT_FILE',
                         help='Optional output name for SVG format report file. This option will be ignored if more than two paired FASTA files are input (i.e. for batch mode); automated naming will be used instead.')

  arg_parse.add_argument('-b', metavar='EXE_FILE',
                         help='Path to bowtie2 (read aligner) executable (will be searched for if not specified)')

  arg_parse.add_argument('-q', metavar='SCHEME',
                         help='Use a specific FASTQ quality scheme (normally not set and deduced automatically). ' + avail_quals)

  arg_parse.add_argument('-qm', default=DEFAULT_MIN_QUAL, metavar='MIN_QUALITY', type=int,
                         help='Minimum acceptable FASTQ quality score in range 0-40 for clipping 3\' end of reads. Default: %d' % DEFAULT_MIN_QUAL)

  arg_parse.add_argument('-m', default=False, action='store_true',
                         help='Force a re-mapping of genome restriction enzyme sites (otherwise cached values will be used if present)')

  arg_parse.add_argument('-p', default=False, action='store_true',
                         help='The input data is multi-cell/population Hi-C; single-cell processing steps are avoided')

  arg_parse.add_argument('-pt', nargs=2, metavar='PAIRED_READ_TAGS', default=['r_1','r_2'],
                         help='When more than two FASTQ files are input (batch mode), the subtrings/tags which differ between paired FASTQ file paths. Default: r_1 r_2')

  arg_parse.add_argument('-x', '--reindex', default=False, action='store_true', dest='x',
                         help='Force a re-indexing of the genome (given appropriate FASTA files)')

  arg_parse.add_argument('-f', nargs='+', metavar='FASTA_FILES',
                         help='Specify genome FASTA files for genome index building (accepts wildcards)')

  arg_parse.add_argument('-f2', nargs='+', metavar='FASTA_FILES_2',
                         help='A second set of genome FASTA files for building a second genome index when using hybrid strain cells (accepts wildcards).')

  arg_parse.add_argument('-a', default=False, action='store_true',
                         help='Whether to report ambiguously mapped contacts')

  arg_parse.add_argument('-k', default=False, action='store_true',
                         help='Keep any intermediate files (e.g. clipped FASTQ etc).')

  arg_parse.add_argument('-sam', default=False, action='store_true',
                         help='Write paired contacts files to SAM format')

  arg_parse.add_argument('-l', metavar='SEQUENCE',
                         help='Seek a specific ligation junction sequence (otherwise this is guessed from the primary restriction enzyme)')

  arg_parse.add_argument('-z', default=False, action='store_true',
                         help='GZIP compress any output FASTQ files')

  arg_parse.add_argument('-v', '--verbose', action='store_true', dest='v',
                         help='Display verbose messages to report progress')

  arg_parse.add_argument('-hc', '--homologous_chromos', dest='hc', metavar='HOM_CHROMO_TSV_FILE',
                         help='File path specifying whitespace-separated pairs of homologous chromosome names (to match first word in header lines of genome sequence FASTQ files) and corresponding pair name (e.g. "chrX"). Required forhybrid strain analysis. See -g2 option.')

  arg_parse.add_argument('-u', default=False, action='store_true',
                         help='Whether to only accept uniquely mapping genome positions and not attempt to resolve certain classes of ambiguous mapping where a single perfect match is found,')

  arg_parse.add_argument('-c', default=0, metavar='GENOME_COPIES',
                         type=int, help='Number of whole-genome copies, e.g. for S2 phase; Default 1 unless homologous chromosomes (-hc) are specified for hybrid genomes.')

  args = vars(arg_parse.parse_args(argv))

  fastq_paths = args['i']
  genome_index = args['g']
  genome_index2 = args['g2']
  re1 = args['re1']
  re2 = args['re2']
  sizes = args['s']
  num_cpu = args['n']
  min_rep = args['r']
  ambig = args['a']
  out_file = args['o']
  ambig_file = args['oa']
  report_file = args['or']
  align_exe = args['b']
  qual_scheme = args['q']
  min_qual = args['qm']
  g_fastas = args['f']
  g_fastas2 = args['f2']
  is_pop_data = args['p']
  pair_tags = args['pt']
  remap = args['m']
  reindex = args['x']
  keep_files = args['k']
  lig_junc = args['l']
  zip_files = args['z']
  unique_map = args['u']
  sam_format = args['sam']
  verbose = args['v']
  homo_chromo = args['hc']
  num_copies = args['c']

  if not num_copies:
    if homo_chromo:
      num_copies = 2
    else:
      num_copies = 1

  if sizes:
    sizes = sorted([int(x) for x in re.split('\D+', sizes)])

  if not fastq_paths:
    msg = 'No FASTQ paths specified'
    fatal(msg)

  elif len(fastq_paths) > 2: # Batch mode
    fastq_paths_1, fastq_paths_2 = pair_fastq_files(fastq_paths, pair_tags)

    if out_file or ambig_file or report_file:
      msg = 'Output naming options ignored for batch mode. Automated naming will be used for each FASTQ pair'
      warn(msg)

    for fastq_path_pair in zip(fastq_paths_1, fastq_paths_2):
      nuc_process(fastq_path_pair, genome_index, genome_index2, re1, re2, sizes, min_rep,
                  num_cpu, num_copies, ambig, unique_map, homo_chromo, None, None, None, align_exe,
                  qual_scheme, min_qual, g_fastas, g_fastas2, is_pop_data, remap, reindex, keep_files,
                  lig_junc, zip_files, sam_format, verbose)

  else:
    nuc_process(fastq_paths, genome_index, genome_index2, re1, re2, sizes, min_rep, num_cpu, num_copies,
                ambig, unique_map, homo_chromo, out_file, ambig_file, report_file, align_exe,
                qual_scheme, min_qual, g_fastas, g_fastas2, is_pop_data, remap, reindex, keep_files,
                lig_junc, zip_files, sam_format, verbose)

  # Required:
  #  - Output CSV report file option
  #
  # Next
  #  - Documentation
  #    + Main PDF/HTML, brief tutorial example
  #
  # To think about:
  #  - could add an option to separate isolated contacts - needs Python only version
  #  - options for different RE strategies, e.g. no fill-in of sticky ends etc.
  #  - normalise/correct population data?
  #  - split input files to parallelise the initial stages of the process


if __name__ == '__main__':
  main()
