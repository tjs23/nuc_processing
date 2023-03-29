"""
---- COPYRIGHT ----------------------------------------------------------------

Copyright (C) 20016-2022
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
from io import BufferedReader

from collections import defaultdict
from shutil import move
from subprocess import Popen, PIPE, call
from math import floor

PROG_NAME = 'nuc_process'
VERSION = '1.3.2'
DESCRIPTION = 'Chromatin contact paired-read Hi-C processing module for Nuc3D and NucTools'
RE_CONF_FILE = 'enzymes.conf'
RE_SITES = {'MboI'   : '^GATC_',
            'DpnII'  : '^GATC_',
            'DpnII*' : '^GATC_',
            'AluI'   : 'AG^CT',
            'BglII'  : 'A^GATC_T',
            'HindIII': 'A^AGCT_T',
            'HindIII*': 'A^AGCT_T'}
STAR_SITES = {'HindIII*':('A^AACT_T',
                          'A^ACCT_T',
                          'A^AGAT_T',
                          'A^AGCA_T',
                          'A^AGCC_T',
                          'A^AGCG_T',
                          'A^AGCT_A',
                          'A^AGCT_C',
                          'A^AGCT_G'),
              'DpnII*':('^GATG_',
                        '^CATG_',
                        '^CATC_')}
ADAPTER_SEQS = {'Nextera':'CTGTCTCTTATA',
                'Illumina universal':'AGATCGGAAGAGC'}
DEFAULT_ADAPTER = 'Illumina universal'
QUAL_SCHEMES = ['phred33', 'phred64', 'solexa']
DEFAULT_MIN_QUAL = 10
QUAL_ZERO_ORDS = {'phred33':33, 'phred64':64, 'solexa':64}
FASTQ_READ_CHUNK = 1048576
READ_BUFFER = 2**16
MIN_READ_LEN = 18
NUM_MAP_FASTAS = 10
CLOSE_AMBIG = 1000
BOWTIE_MAX_AMBIG_SCORE_TOL = 5
EXCLUDE_BIN_SIZE = 100000
ID_LEN = 10
MIN_ADAPT_OVERLAP = 7
INF = float('inf')
SCORE_TAG = re.compile(r'\sAS:i:(\S+).+\sMD:Z:(\S+)')
SCORE_TAG_SEARCH = SCORE_TAG.search
FILENAME_SPLIT_PATT = re.compile('[_\.]')
NCC_FORMAT = '%s %d %d %d %d %s %s %d %d %d %d %s %.1f %d %d\n'
LOG_FILE_PATH = None
STAT_FILE_PATH = None
STAT_SECTIONS = {'command':'Command',
                 'general':'General',
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


def open_file_r(file_path, complete=True, gzip_exts=('.gz','.gzip'), buffer_size=READ_BUFFER):

  import io, subprocess
  
  if os.path.splitext(file_path)[1].lower() in gzip_exts:

    if complete:
      try:
        file_obj = subprocess.Popen(['zcat', file_path], stdout=subprocess.PIPE).stdout
      except OSError:
        file_obj = BufferedReader(gzip.open(file_path, 'rb'), buffer_size)
    
    else:
      file_obj = BufferedReader(gzip.open(file_path, 'rb'), buffer_size)
    
    if sys.version_info.major > 2:
      file_obj = io.TextIOWrapper(file_obj, encoding="utf-8")
 
  else:
    if sys.version_info.major > 2:
      file_obj = open(file_path, 'rU', buffer_size, encoding='utf-8')
      
    else:
      file_obj = open(file_path, 'rU', buffer_size)
  
  return file_obj
 


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

     id1 = int(line1[:ID_LEN])
     id2 = id1-1

     while (id2 < id1) and line2:
       id2 = int(line2[:ID_LEN])
       line2 = read2()

     if (id1 == id2) and (id1 in sam_ids):
       if sam_ids[id1]:
         line2, line1 = line1, line2

       write('%s\n%s\n' % pair_sam_lines(line1[ID_LEN:], line2[ID_LEN:]))
   
   in_file_obj_1.close()
   in_file_obj_2.close()
   sam_file_obj.close()
   ncc_file_obj.close()

   return sam_file_path


def merge_file_names(file_path1, file_path2, sep='_'):

  # same dir, need non truncated name

  dir_name1, file_name1 = os.path.split(file_path1)
  dir_name2, file_name2 = os.path.split(file_path2)

  if dir_name1 != dir_name2:
    msg = 'Attempt to merge file names for file from different directories'
    raise Exception(msg)

  file_root1, file_ext1 = os.path.splitext(file_name1)
  file_root2, file_ext2 = os.path.splitext(file_name2)

  if file_ext1 != file_ext2:
    msg = 'Attempt to merge file names with different file extensions'
    raise Exception(msg)

  parts1 = FILENAME_SPLIT_PATT.split(file_root1)
  parts2 = FILENAME_SPLIT_PATT.split(file_root2)
  parts3 = []

  n1 = len(parts1)
  n2 = len(parts2)
  n = max(n1, n2)

  for i in range(n):

    if (i < n1) and (i < n2):
      a = parts1[i]
      b = parts2[i]

      parts3.append(a)
      if a != b:
        parts3.append(b)

    elif i < n1:
      parts3.append(parts1[i])
    else:
      parts3.append(parts2[i])

  file_root3 = sep.join(parts3)

  file_path3 = os.path.join(dir_name1, file_root3 + file_ext1)

  return file_path3


def remove_promiscuous(ncc_file, num_copies=1, keep_files=True, zip_files=False,
                       resolve_limit=1e3, close_cis=1e4, ambig=False):
  """
  Function to remove contacts with promiscuous ends from an NCC format file.
  Promiscuous ends occur where a specififc restriction fragment end is involved
  in more contacts that would normally be allowed given the ploidy of the cell.

  resolve_limit : Allow two suitably close promiscous ends if the pairs are long range cis or trans
  """
  
  from itertools import combinations
  
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
  for line in in_file_obj:
    line_data = line.split()
    chr_a = line_data[0]
    chr_b = line_data[6]
    f_start_a = int(line_data[3])
    f_start_b = int(line_data[9])
    strand_a = line_data[5]
    strand_b = line_data[11]
    ambig_sz = int(float(line_data[12]))
    read_group = int(line_data[13])
    frag_counts[(chr_a, f_start_a, strand_a)].add((chr_b, f_start_b, strand_b, read_group))
    frag_counts[(chr_b, f_start_b, strand_b)].add((chr_a, f_start_a, strand_a, read_group))

  remove_read_group = set()
  resolved_read_group = set()
  
  if hasattr(frag_counts, 'iteritems'):
    items = frag_counts.iteritems
  else:
    items = frag_counts.items

  for re_start, re_ends in items():
    if len(re_ends) > num_copies:
      read_groups = {x[3] for x in re_ends}

      if len(read_groups) <= num_copies: # Multiple ends have mo more ambiguity groups than there are genomes/haplotypes
        continue

      if len(re_ends) == num_copies+1:
        chr_a, f_start_a, strand_a = re_start
        
        for re_end1, re_end2 in combinations(re_ends, 2):
          chr_b1, f_start_b1, strand_b1, read_b1 = re_end1
          chr_b2, f_start_b2, strand_b2, read_b2 = re_end2

          if (chr_b1 == chr_b2) and abs(f_start_b1-f_start_b2) < resolve_limit: # Other ends are very close to each other
            if chr_a != chr_b1: # Trans to this end
              resolved_read_group.update(read_groups) # Two groups are effectively the same at one end

            elif abs(f_start_a-f_start_b1) > close_cis and abs(f_start_a-f_start_b2) > close_cis: # Far from this end
              resolved_read_group.update(read_groups) # Two groups are effectively the same at one end

            else: # Other end too close to this end
              remove_read_group.update(read_groups)
              
          else: # Other ends too far apart
            remove_read_group.update(read_groups)

      else: # Cannot resolve with more than two ends
        remove_read_group.update(read_groups) # All associated ambiguity groups annuled

  n_promiscuous = 0
  n_resolved = 0
  n_clean = 0

  in_file_obj.seek(0)
  for line in in_file_obj:
    line_data = line.split()
    chr_a = line_data[0]
    chr_b = line_data[6]
    f_start_a = int(line_data[3])
    f_start_b = int(line_data[9])
    strand_a = line_data[5]
    strand_b = line_data[11]
    ambig_sz = int(float(line_data[12]))
    read_group = int(line_data[13])

    if read_group in remove_read_group:
      if ambig_sz:
        n_promiscuous += 1

      if keep_files:
        write_promisc(line)
    
    elif read_group in resolved_read_group:
      if ambig_sz:
        n_resolved += 1
     
      write_clean(line)
     
    else:
      if ambig_sz:
        n_clean += 1
        
      write_clean(line)

  in_file_obj.close()

  if keep_files:
    if zip_files:
      compress_file(promiscous_ncc_file)

  n = n_promiscuous + n_clean + n_resolved

  stats = [('input_contacts', n),
           ('clean',(n_clean, n)),
           ('promiscuous',(n_promiscuous, n)),
           ('resolved',(n_resolved, n)),
           ('accepted',(n_clean+n_resolved, n)),
           ]

  stat_key = 'promsic_ambig' if ambig else 'promsic'
  log_report(stat_key, stats)

  move(clean_ncc_file_temp, clean_ncc_file)

  return clean_ncc_file


def get_ncc_stats(ncc_file, hom_chromo_dict, far_min=10000):

  n_pairs = 0
  n_ambig_pairs = 0
  n_ambig = 0
  n_contacts = 0
  group_counts = [0,0,0,0] # near, far, trans, homolog
  
  group = set()
  group_add = group.add
  
  with open(ncc_file) as in_file_obj:
    for line in in_file_obj:
      line_data = line.split()
      chr_a = line_data[0]
      chr_b = line_data[6]
      ambig, weight = line_data[12].split('.')
      ambig = int(ambig)
      
      n_pairs += 1
      
      if ambig > 0:
        n_contacts += 1
        
        if ambig > 1:
          n_ambig += 1
          n_ambig_pairs += 1
          
        # Either group is unary or the code is prioritised for homolog > trans > far > near
        # - cold be aonly a single code index for a consistent ambigous group
        if group:
          group_counts[max(group)] += 1
          group = set()
          group_add = group.add
      
      else:
        n_ambig_pairs += 1
      
      if chr_a != chr_b:
        if chr_b == hom_chromo_dict.get(chr_a):
          group_add(3)
        else:
          group_add(2)

      else:
        idx1 = 4 if line_data[5] == '+' else 3  # + strand : Ligation junction is ahead at 3' end of RE frag
        idx2 = 10 if line_data[11] == '+' else 9

        p1 = int(line_data[idx1])
        p2 = int(line_data[idx2])

        if abs(p1-p2) < far_min:
          group_add(0)
        else:
          group_add(1)
  
  if group:
    group_counts[max(group)] += 1

  stats = [('total_contacts', n_contacts),
           ('ambig_contacts', (n_ambig, n_contacts)),
           ('total_pairs', n_pairs),
           ('mean_ambiguity', n_ambig_pairs/float(n_ambig or 1)),
           ('cis_near',(group_counts[0], n_contacts)),
           ('cis_far',(group_counts[1], n_contacts)),
           ('trans',(group_counts[2], n_contacts)),
           ('homolog_trans',(group_counts[3], n_contacts)),
           ]

  return stats


def remove_redundancy(ncc_file, keep_files=True, zip_files=False, min_repeats=2,
                      use_re_fragments=True, is_hybrid=False):

  """
  A:B B:A redundancey taken care of at this stage because the NCC file pairs are internally sorted

  Option to remove redundancy at he RE fragment level (or otherwise at the read level)

  Choose the longest read (not most common?) as the representitive for a merge
  """

  out_file_name = tag_file_name(ncc_file, 'multi_read')
  sort_file_name = tag_file_name(ncc_file, 'sort')

  if INTERRUPTED and os.path.exists(out_file_name) and not os.path.exists(sort_file_name):
    return out_file_name

  # Calculate sizes of ambiguity groups

  ambig_sizes = {}
  n_contacts = 0
  for line in open_file_r(ncc_file):
    ambig, read_id = line.split()[12:14]
    ambig = int(float(ambig,))
    
    if ambig > 0:
      ambig_sizes[read_id] = ambig
      n_contacts += 1

  # Make temporary sorted file
  cmd_args = ['sort', ncc_file]
  call(cmd_args, shell=False, stderr=None, stdin=None, stdout=open(sort_file_name, 'w'))
  sort_file_obj = open(sort_file_name, 'r')

  n_unique = 0
  n_redundant = 0
  mean_redundancy = 0.0

  # Compare using strand information given we want to know which fragment END is used
  if use_re_fragments:
    ncc_idx = (0,3,5,6,9,11) # Compare on RE fragment: chr_a, f_start_a, strand_a, chr_b, f_start_b, strand_b
  else:
    ncc_idx = (0,1,5,6,7,11) # Compare on read starts - does not consider ends as these can vary in read replicates

  # Remove repeats

  excluded_reads = set()
  group_reps = defaultdict(int)
  sort_file_obj.seek(0)
  line_keep = sort_file_obj.readline()

  if line_keep: # Could be empty
    line_data = line_keep.split()

    n = 1
    keep_data = [line_data[i] for i in ncc_idx]
    len_keep = abs(int(line_data[2]) - int(line_data[2])) + abs(int(line_data[8]) - int(line_data[7]))
    read_keep = line_data[13]
    read_ambig_keep = ambig_sizes[read_keep]

    # Normally remove redundancy, keeping the least ambiguous or else the longest
    # - for hybrid dual-genome mapping go for the most ambiguous to be safe
    for line in sort_file_obj:
      line_data = line.split()
      curr_data = [line_data[i] for i in ncc_idx]
      len_curr = abs(int(line_data[2]) - int(line_data[1])) + abs(int(line_data[8]) - int(line_data[7]))
      read_curr = line_data[13]
      read_ambig_curr = ambig_sizes[read_curr]

      if curr_data == keep_data:
        n += 1

        if read_ambig_curr == read_ambig_keep: # Equally ambiguous
          if len_curr > len_keep: # This, longer repeat is better
            excluded_reads.add(read_keep) # Remove prev kept group
            len_keep = len_curr
            read_ambig_keep = read_ambig_curr
            read_keep = read_curr

          else: # This repeat is worse
            excluded_reads.add(read_curr) # Remove this group

        elif read_ambig_curr > read_ambig_keep and is_hybrid: # For hybrids keep the more ambiguous
          excluded_reads.add(read_keep) # Remove prev kept group
          len_keep = len_curr
          read_ambig_keep = read_ambig_curr
          read_keep = read_curr

        elif read_ambig_curr < read_ambig_keep and not is_hybrid: # Normally better to keep this less ambiguous repeat
          excluded_reads.add(read_keep) # Remove prev kept group
          len_keep = len_curr
          read_ambig_keep = read_ambig_curr
          read_keep = read_curr

        else: # Keep the old one
          excluded_reads.add(read_curr) # Remove this group

        continue
     
      else: # Not a repeat, write previous
        group_reps[read_keep] += n
        mean_redundancy += n
        len_keep = len_curr
        keep_data = curr_data
        read_keep = read_curr
        read_ambig_keep = read_ambig_curr
        n = 1

  if keep_files:
    uniq_file_name = tag_file_name(ncc_file, 'unique_read')
    uniq_file_obj = open(uniq_file_name, 'w')
    uniq_write = uniq_file_obj.write
  
  # Remove excluded ambig groups and non-supported, preserving original file order

  n_excluded = len(excluded_reads)
  
  with open(ncc_file) as in_file_obj, open(out_file_name, 'w') as out_file_obj:
    write = out_file_obj.write
    
    if use_re_fragments: # E.g. single-cell
      for line in in_file_obj:
        ambig, read_id = line.split()[12:14]
        ambig = int(float(ambig))

        if read_id in excluded_reads:
          continue

        elif group_reps[read_id]/float(ambig_sizes[read_id]) < min_repeats:
          if keep_files:
            uniq_write(line)
          
          if ambig > 0: 
            n_unique += 1

        else:
          write(line)
          
          if ambig > 0: 
            n_redundant += 1
    
    else: # E.g. bulk/population
      for line in in_file_obj:
        ambig, read_id = line.split()[12:14]
        ambig = int(float(ambig))

        if read_id in excluded_reads:
          continue
        
        if ambig > 0:
          if group_reps[read_id]/float(ambig_sizes[read_id]) < 2:
             n_unique += 1
          else:
             n_redundant += 1
          
        write(line)
  
  if keep_files:
    if zip_files:
      compress_file(ncc_file)

  else:
    os.unlink(ncc_file)

  os.unlink(sort_file_name) # Remove temp sorted file
  
  if keep_files:
    uniq_file_obj.close()

  n = n_unique + n_redundant
  mean_redundancy /= float(n) or 1.0

  stats = [('input_contacts', n_contacts),
           ('redundant_contacts', (n_excluded, n_contacts)),
           ('effective_contacts', n),
           ('unique', (n_unique, n)),
           ('supported', (n_redundant, n)),
           ('mean_redundancy', mean_redundancy)]
  
  log_report('dup', stats)

  return out_file_name


def filter_pairs(pair_ncc_file, re1_files, re2_files, chromo_name_dict, hom_chromo_dict,
                 sizes=(100,2000), keep_files=True, zip_files=False, min_mappability=1,
                 re2_tolerance=1000, ambig=False, star_dict=None):

  filter_file = tag_file_name(pair_ncc_file, 'filter_accepted', '.ncc')
  filter_file_temp = filter_file + TEMP_EXT

  if INTERRUPTED and os.path.exists(filter_file) and not os.path.exists(filter_file_temp):
    return filter_file, {}

  re1_frag_dict = {}
  re1_end_dict = {}
  re2_frag_dict = {}
  re2_end_dict = {}

  for i, re1_file in enumerate(re1_files):
    frags, ends = read_re_frag_file(re1_file, i, hom_chromo_dict)
    re1_frag_dict.update(frags)
    re1_end_dict.update(ends)

  for i, re2_file in enumerate(re2_files):
    frags, ends = read_re_frag_file(re2_file, i, hom_chromo_dict)
    re2_frag_dict.update(frags)
    re2_end_dict.update(ends)
    
  out_file_objs = {}
  write_funcs = {}
  out_file_names = {}
  counts = {}
  min_size, max_size = sorted(sizes)
  size_counts = defaultdict(int)
  size_counts_accept_cis = defaultdict(int)
  size_counts_accept_trans = defaultdict(int)
  junct_sep_counts_pos = defaultdict(int)
  junct_sep_counts_neg = defaultdict(int)
  
  valid_chromos = set(chromo_name_dict.values())

  if re1_files:
    valid_chromos &= set(re1_end_dict.keys())
    use_re = True
    too_close_limit = int((min_size+max_size)/2)
  else:
    use_re = False
    too_close_limit = 75 # less than half a nucleosome
  
  
  for tag in ('accepted', 'too_small', 'too_big', 'circular_re1', 'internal_re1',
              'internal_re2', 'no_end_re2', 'overhang_re1', 'too_close',
              'adjacent_re1', 'unknown_contig', 'excluded_group', 'no_insert'):
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
  
  def adjust_cut_site(pos_strand, read_start, re_start, re_end, star_pos):

     if pos_strand:
       re_end = star_pos
       delta_re = max(0, re_end - read_start)
       
     else:
       re_start = star_pos
       delta_re = max(0, read_start - re_start)
    
     return re_start, re_end, delta_re
    
  in_file_obj = open_file_r(pair_ncc_file)
  excluded_reads = set()
  pruned_groups = set()
  
  if re1_files:
    max_idx = {chr_a:len(re1_end_dict[chr_a])-1 for chr_a in re1_end_dict}
  
  if re2_files:
    max_idx2 = {chr_a:len(re2_end_dict[chr_a])-1 for chr_a in re2_end_dict}
  
  if star_dict:
    max_star_idx = {chr_a:len(star_dict[chr_a])-1 for chr_a in star_dict}
  
  searchsorted = np.searchsorted
  
  for line in in_file_obj:
    # NCC file starts are always the start of the READ, not what was in the BAM file
    # which means start > end where strand is "-"

    chr_a, start_a, end_a, f_start_a, f_end_a, strand_a, chr_b, start_b, end_b, \
      f_start_b, f_end_b, strand_b, ambig_code, read_id, swap_pair = line.split()

    counts['input_pairs'] += 1

    pos_strand_a = strand_a == '+'
    pos_strand_b = strand_b == '+'
    
    start_a = int(start_a)
    start_b = int(start_b)
    
    end_a = int(end_a)
    end_b = int(end_b)
    
    read_id = int(read_id)
    
    if chr_a not in valid_chromos:
      count_write('unknown_contig', line)
      excluded_reads.add(read_id)
      continue

    if chr_b not in valid_chromos:
      count_write('unknown_contig', line)
      excluded_reads.add(read_id)
      continue

    if re2_files:
      if chr_a not in re2_end_dict:
        count_write('unknown_contig', line)
        excluded_reads.add(read_id)
        continue

      if chr_b not in re2_end_dict:
        count_write('unknown_contig', line)
        excluded_reads.add(read_id)
        continue
    
    is_cis = chr_a == chr_b

    # Find which RE fragments the _final_ (3') read positions were within ; this points to the following RE1 ligation junction
    # - the read could cover undigested RE sites
    # - index of frag with pos immediately less than end
    
    # With truncated contigs the (real) read sequence can exceed the reference
    if use_re:
      re1_a_idx = min(max_idx[chr_a], searchsorted(re1_end_dict[chr_a], end_a))
      re1_b_idx = min(max_idx[chr_b], searchsorted(re1_end_dict[chr_b], end_b))
 
      re1_a_start, re1_a_end, mappability_a = re1_frag_dict[chr_a][re1_a_idx]
      re1_b_start, re1_b_end, mappability_b = re1_frag_dict[chr_b][re1_b_idx]

      if pos_strand_a:
        delta_re1_a = max(0, re1_a_end - start_a) # separation from ligation junction
        p1, p2 = start_a, end_a # sorted GENOME positions
      else:
        delta_re1_a = max(0, start_a - re1_a_start)
        p1, p2 = end_a, start_a

      if pos_strand_b:
        delta_re1_b = max(0, re1_b_end - start_b)
        p3, p4 = start_b, end_b
      else:
        delta_re1_b = max(0, start_b - re1_b_start)
        p3, p4 = end_b, start_b

      size_t = delta_re1_a + delta_re1_b

      # With some REs (e.g. HindIII) fragments that are apprently too big may be due to star activity
 
      if star_dict and size_t > max_size and len(star_dict[chr_a]) and len(star_dict[chr_b]):

        star_a_idx = searchsorted(star_dict[chr_a], start_a)
        star_b_idx = searchsorted(star_dict[chr_b], start_b)
 
        if not pos_strand_a:
          star_a_idx -= 1

        if not pos_strand_b:
          star_b_idx -= 1
 
        star_a_idx = max(0, min(max_star_idx[chr_a], star_a_idx))
        star_b_idx = max(0, min(max_star_idx[chr_b], star_b_idx))
 
        star_a_pos = star_dict[chr_a][star_a_idx]
        star_b_pos = star_dict[chr_b][star_b_idx]
 
        if re1_a_start < star_a_pos < re1_a_end:
          if re1_b_start < star_b_pos < re1_b_end: # Potential star activity both sides
            if delta_re1_a > delta_re1_b: # Trim longest side first
              re1_a_start, re1_a_end, delta_re1_a = adjust_cut_site(pos_strand_a, start_a, re1_a_start, re1_a_end, star_a_pos)
 
              if delta_re1_a + delta_re1_b > max_size:
                re1_b_start, re1_b_end, delta_re1_b = adjust_cut_site(pos_strand_b, start_b, re1_b_start, re1_b_end, star_b_pos)
 
            else:
              re1_b_start, re1_b_end, delta_re1_b = adjust_cut_site(pos_strand_b, start_b, re1_b_start, re1_b_end, star_b_pos)
 
              if delta_re1_a + delta_re1_b > max_size:
                re1_a_start, re1_a_end, delta_re1_a = adjust_cut_site(pos_strand_a, start_a, re1_a_start, re1_a_end, star_a_pos)
 
          else:
            re1_a_start, re1_a_end, delta_re1_a = adjust_cut_site(pos_strand_a, start_a, re1_a_start, re1_a_end, star_a_pos)
 
        elif re1_b_start < star_b_pos < re1_b_end:
          re1_b_start, re1_b_end, delta_re1_b = adjust_cut_site(pos_strand_b, start_b, re1_b_start, re1_b_end, star_b_pos)
 
        size_t = delta_re1_a + delta_re1_b # Re-check

    else: # No RE digestion, e.g. MicroC
    
      if pos_strand_a:
        p1, p2 = start_a, end_a # sorted GENOME positions
      else:
        p1, p2 = end_a, start_a

      if pos_strand_b:
        p3, p4 = start_b, end_b
      else:
        p3, p4 = end_b, start_b
      
      size_t = min_size + (end_a - start_a) + (end_b - start_b)
      re1_a_start, re1_a_end = start_a, end_a
      re1_b_start, re1_b_end = start_b, end_b
      re1_a_idx, re1_b_idx = 0, 2 # avoids same-RE checks
      
      if is_cis:
        delta_re1_a =  delta_re1_b = abs(end_a-end_b)
      else:
        delta_re1_a = 0
        delta_re1_b = 0
      
    # Collect stats
    size_ok = min_size <= size_t <= max_size
    

    if use_re:
      size_bin = int(10*np.log10(max(size_t, 1)))
      
      if pos_strand_a:
        junct_sep_counts_pos[int(delta_re1_a/10.0)] += 1
      else:
        junct_sep_counts_neg[int(delta_re1_a/10.0)] += 1
        
      if pos_strand_b:
        junct_sep_counts_pos[int(delta_re1_b/10.0)] += 1
      else:
        junct_sep_counts_neg[int(delta_re1_b/10.0)] += 1
    
    elif is_cis:
      size_bin = int(10*np.log10(max(delta_re1_a, 1)))
      
      if strand_a == strand_b: # Same strand
        junct_sep_counts_pos[int(delta_re1_a/10.0)] += 1
      else:
        junct_sep_counts_neg[int(delta_re1_a/10.0)] += 1
      
    else:
      size_bin = int(10*np.log10(max(size_t, 1)))
    
    size_counts[size_bin] += 1
    
    # Add RE fragment positions
    line = NCC_FORMAT % (chr_a, start_a, end_a, re1_a_start, re1_a_end, strand_a,
                         chr_b, start_b, end_b, re1_b_start, re1_b_end, strand_b,
                         float(ambig_code), read_id, int(swap_pair))
    
    if re2_files:
      # Check Read is at the end of different RE2 fragments
      # Note: modulo is due to circular chromosomes
      re2_a_idx = min(max_idx2[chr_a], searchsorted(re2_end_dict[chr_a], end_a))
      re2_b_idx = min(max_idx2[chr_b], searchsorted(re2_end_dict[chr_b], end_b))

      # Internal RE2 not necessarily a problem if reads are in different RE1 fragements, e.g. potential ligation junction within
      if re2_a_idx == re2_b_idx and chr_a == chr_b and re1_a_idx == re1_b_idx:
        count_write('internal_re2', line)
        excluded_reads.add(read_id) # The real pair could definitely be a duff one
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
        excluded_reads.add(read_id) # The real pair could definitely be a duff one
        #pruned_groups.add(read_id)
        continue

      if abs(delta_re2_b) > re2_tolerance:
        count_write('no_end_re2', line)
        excluded_reads.add(read_id) # The real pair could definitely be a duff one
        #pruned_groups.add(read_id)
        continue

    if is_cis:
      if re1_a_idx == re1_b_idx: # Same Re1 fragment
        if start_a < start_b:
          if pos_strand_b and not pos_strand_a: # Reads go outwards in same frag: <-A B->
            count_write('circular_re1', line)
            excluded_reads.add(read_id)
            continue

        else:
          if pos_strand_a and not pos_strand_b: # Reads go outwards in same frag: <-B A->
            count_write('circular_re1', line)
            excluded_reads.add(read_id)
            continue

        if (p1 < re1_a_start < p2) or (p1 < re1_a_end < p2): # Read A overshoots Re1 fragment end
          count_write('overhang_re1', line)

        elif (p3 < re1_b_start < p4) or (p3 < re1_b_end < p4): # Read B overshoots Re1 fragment end
          count_write('overhang_re1', line)

        else:
          count_write('internal_re1', line)

        excluded_reads.add(read_id) # The real pair could definitely be a duff one

      elif abs(re1_a_idx-re1_b_idx) < 2:
        count_write('adjacent_re1', line) # Mostly re-ligation
        excluded_reads.add(read_id) # The real pair could definitely be useless

      else: # Different Re1 fragment or no Re1

        if pos_strand_a != pos_strand_b:

          if pos_strand_a and (p1 < p3): # Sequencing toward each other
            delta = p4 - p1 # separation of pair
            
            if delta < too_close_limit: # Pair is possible even without ligation
              if delta < 1:
                count_write('no_insert', line) 
              else:
                count_write('too_close', line) 

              excluded_reads.add(read_id)  # The real pair could definitely be useless
              continue

          elif pos_strand_b and (p3 < p1): # Sequencing toward each other
            delta = p2 - p3

            if delta < too_close_limit: # Pair is possible even without ligation
              if delta < 1:
                count_write('no_insert', line) 
              else:
                count_write('too_close', line) 
                
              excluded_reads.add(read_id) # The real pair could definitely be useless
              continue

        if size_ok:
          if min(abs(p1-p3), abs(p2-p4)) < 1e4:
            counts['near_cis_pairs'] += 1
          else:
            counts['far_cis_pairs'] += 1

          count_write('accepted', line)
          size_counts_accept_cis[size_bin] += 1

        elif size_t < min_size:
          count_write('too_small', line)
          excluded_reads.add(read_id)

        else:
          count_write('too_big', line)
          excluded_reads.add(read_id)

    else: # Good stuff, trans contact, many errors impossible
      if size_ok:
        counts['trans_pairs'] += 1
        count_write('accepted', line)
        size_counts_accept_trans[size_bin] += 1

      elif size_t < min_size:
        count_write('too_small', line)
        excluded_reads.add(read_id)

      else:
        count_write('too_big', line)
        
        #if size_t < 2 * max_size: # Not impossible
        excluded_reads.add(read_id)

  # Remove complete excluded ambiguity groups at the end:
  # - Removes a whole ambiguity group if only one possibilty was suspect as it could be the real contact
  
  in_file_obj.close()
  out_file_objs['accepted'].close()  
  del out_file_objs['accepted']
  prune_counts = defaultdict(int)
  
  if pruned_groups: # Degree of ambiguity changed
    with open_file_r(out_file_names['accepted']) as file_obj:
      for line in file_obj:
        read_id = int(line.split()[13])
        
        if read_id in pruned_groups:
          prune_counts[read_id] += 1
  
  out_file_obj = open(filter_file, 'w')
  with open_file_r(out_file_names['accepted']) as file_obj:
    write = out_file_obj.write

    for line in file_obj:
      read_id = int(line.split()[13])

      if read_id in excluded_reads:
        count_write('excluded_group', line)

      else:
        if read_id in prune_counts:
          count = prune_counts[read_id]
          row = line.split()
          row[12] = '%.1f' % (count + 0.1)
          line = ' '.join(row) + '\n'
          prune_counts[read_id] = 0.0
        
        write(line)

  counts['accepted'] -= counts['excluded_group']

  for tag in out_file_objs:
    out_file_objs[tag].close()

  n = counts['input_pairs']
  stats_list = []
  
  if use_re:
    stat_keys = ('input_pairs', 'accepted', 'near_cis_pairs', 'far_cis_pairs', 'trans_pairs',
                 'internal_re1', 'adjacent_re1', 'circular_re1', 'overhang_re1', 'no_insert',
                 'too_close', 'too_small', 'too_big', 'excluded_group',
                 'internal_re2', 'no_end_re2', 'unknown_contig')
    
  else:
    stat_keys = ('input_pairs', 'accepted', 'near_cis_pairs', 'far_cis_pairs', 'trans_pairs',
                 'no_insert', 'too_close', 'too_small', 'too_big', 'excluded_group', 'unknown_contig')
  
  for key in stat_keys:
    stats_list.append((key, (counts[key], n))) # So percentages can be calculated in reports

  del out_file_names['accepted'] # Main output never compressed at this stage

  if keep_files and zip_files:
    for tag in out_file_objs:
      compress_file(out_file_names[tag])
  
  m = max([max(size_counts_accept_trans.keys() or [0]),
           max(size_counts_accept_cis.keys() or [0]),
           max(size_counts.keys() or [0])])
  
  size_edges = [x/10.0 for x in range(m)]

  n_cis = 0.01 * float(counts['near_cis_pairs'] + counts['far_cis_pairs']) or 1.0
  hist_sizes_accept_cis = [size_counts_accept_cis[i]/n_cis for i in range(m)]
  
  n_trans = 0.01 * float(counts['trans_pairs']) or 1.0
  hist_sizes_accept_trans = [size_counts_accept_trans[i]/n_trans for i in range(m)]
  
  n = n_cis + n_trans
  hist_sizes = [size_counts[i]/n for i in range(m)]

  seps_edges = []  
  histjunct_sep_neg = []
  histjunct_sep_pos = []
  
  for i in range(100):
    seps_edges.append((i+1)*10)
    
    a = junct_sep_counts_pos[i]
    b = junct_sep_counts_neg[i]
    
    if a:
      histjunct_sep_pos.append((a/10.0))  # Was log_10
    else:
      histjunct_sep_pos.append(0.0)
    
    if b:
      histjunct_sep_neg.append((b/10.0))
    else:
      histjunct_sep_neg.append(0.0)
    
  if ambig:
    stat_key = 'filter_ambig'
    frag_key = 'frag_sizes_ambig'
  
  else:
    stat_key = 'filter'
    frag_key = 'frag_sizes'

  extra_vals = (hist_sizes, hist_sizes_accept_cis, hist_sizes_accept_trans, size_edges,
                histjunct_sep_pos, histjunct_sep_neg, seps_edges)

  log_report(stat_key, stats_list, {frag_key:extra_vals})

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


def _write_ncc_line(ncc_data_a, ncc_data_b, ambig_code, read_id, write_func):
  # write NCC line in chromosome, position order for merging etc

  chr_a, start_a = ncc_data_a[:2]
  chr_b, start_b = ncc_data_b[:2]

  if chr_a == chr_b:
    if start_a > start_b:
      line_data = ncc_data_b + ncc_data_a + (ambig_code+0.1, read_id, 1)

    else:
      line_data = ncc_data_a + ncc_data_b + (ambig_code+0.1, read_id, 0)

  elif chr_a > chr_b:
    line_data = ncc_data_b + ncc_data_a + (ambig_code+0.1, read_id, 1)

  else:
    line_data = ncc_data_a + ncc_data_b + (ambig_code+0.1, read_id, 0)

  write_func(NCC_FORMAT % line_data)


def _load_mappability_files(contig_sources, *file_paths):

  reg_dict = defaultdict(set)
  reg_size = None
  
  if not contig_sources:
    contig_sources = {}
  
  for file_path in file_paths:
    with open_file_r(file_path) as file_obj:
      for line in file_obj:
        contig, start, end = line.split()
        source = contig_sources.get(contig, contig)
        start = int(start)

        if not reg_size:
          reg_size = int(end) - start
        
        reg_dict[source].add(start)

  return reg_dict, reg_size


def _load_exclusion_file(file_path, bin_size=EXCLUDE_BIN_SIZE):
  
  exclusion_dict = defaultdict(list)
  
  with open_file_r(file_path) as file_obj:
    for line in file_obj:
      contig, start, end = line.split()
      i = int(start)//bin_size
      j = int(end)//bin_size
      
      key = (contig, i)
      exclusion_dict[key].append((start,end))
      
      if j != i:
        key = (contig, j)
        exclusion_dict[key].append((start,end))
  
  for key in exclusion_dict:
    regions = np.array(sorted(exclusion_dict[key]))
    exclusion_dict[key] = regions[:,0], regions[:,1]
  
  return exclusion_dict
  
       
def pair_mapped_hybrid_seqs(sam_file1, sam_file2, sam_file3, sam_file4, chromo_names,
                            file_root, ambig=True, unique_map=False):

  paired_ncc_file_name = tag_file_name(file_root, 'pair', '.ncc')
  paired_ncc_file_name_temp = paired_ncc_file_name + TEMP_EXT

  if INTERRUPTED and os.path.exists(paired_ncc_file_name) and not os.path.exists(paired_ncc_file_name_temp):
    return paired_ncc_file_name

  ncc_file_obj = open(paired_ncc_file_name_temp, 'w')
  write_pair = ncc_file_obj.write

  file_objs = [open_file_r(sam_file) for sam_file in (sam_file1, sam_file2, sam_file3, sam_file4)]
  readlines = [f.readline for f in file_objs]
      
  n_map = [0,0,0,0]
  n_pairs = 0
  n_ambig = 0
  n_unpaired = 0
  n_unambig = 0
  n_unmapped = 0
  n_hybrid_ambig = 0
  n_hybrid_poor = 0
  n_hybrid_end_missing = 0
  n_primary_strand = [0,0,0,0]
  n_strand = [0,0,0,0]
  
  max_score = 0
  zero_ord = QUAL_ZERO_ORDS['phred33']
  really_bad_score = 2 * (max_score - 2 * (BOWTIE_MAX_AMBIG_SCORE_TOL+1) )
  close_second_best = 2 * BOWTIE_MAX_AMBIG_SCORE_TOL - 1
  
  # Go through same files and pair based on matching id
  # Write out any multi-position mapings

  # Skip to end of headers

  lines = [_skip_sam_header(readline) for readline in readlines]

  # Process data lines

  ids = [line[:ID_LEN] for line in lines]
  searchsorted = np.searchsorted

  while '' not in ids:
    _id = max(ids)
    
    unpaired = set()
    for i in range(4):
      while ids[i] and ids[i] < _id: # Reads may not match as some ends can be removed during clipping
        unpaired.add(ids[i])
        line = readlines[i]()
        lines[i] = line
        ids[i] = line[:ID_LEN]
      
    if unpaired:
      n_unpaired += len(unpaired)
      continue 
      
    else:
      contacts = [[], [], [], []]
      scores = [[], [], [], []]
      
      for i in range(4):
        j = 0
        
        while ids[i] == _id:
          n_map[i] += 1
          line = lines[i]
          data = line.split('\t')
          chromo = data[2]
          chr_name = chromo_names.get(chromo, chromo)
          
          if chromo != '*':
            revcomp = int(data[1]) & 0x10
            start = int(data[3])
            seq = data[9]
            qual = data[10]
            nbp = len(seq)
            end = start + nbp
            
            score, var = SCORE_TAG_SEARCH(line).groups()
            var_orig = var             
            score = int(score) # when --local : -nbp*2

            if revcomp and var[-1] == '0' and var[-2] in 'GCAT': # Ignore substitutions at the end e.g "C0"; add back subtracted score
              var = var[:-2]
            
              if seq[-1] != 'N': # No penalty for N's
                q = min(ord(qual[-1]) - zero_ord, 40.0)
                mp1 = 2 + floor(4*q/40.0)  # MX = 6, MN = 2.
                score = min(score+mp1, max_score)
                end -= 1

            elif var[0] == '0' and var[1] in 'GCAT':
              var = var[2:]
              
              if seq[0] != 'N':
                q = min(ord(qual[0]) - zero_ord, 40.0)
                mp2 = 2 + floor(4*q/40.0)  
                score = min(score+mp2, max_score)
                start += 1
           
            if revcomp:
              start, end = end, start  # The sequencing read started from the other end
            
            ncc = (chr_name, start, end, 0, 0, '-' if revcomp else '+') # Strand info kept because ends can be diffferent for replicate reads, no Re fragment positions, yet
            contacts[i].append((ncc, score))
            scores[i].append(score)
            
            if j == 0:
              n_strand[i] += 1.0
 
              if not revcomp:
                n_primary_strand[i] += 1            
            
          line = readlines[i]()
          lines[i] = line
          ids[i] = line[:ID_LEN]
          j += 1
          
      if not unique_map: # Resolve some perfect vs non-perfect positional ambiguity
        for a in (0,1,2,3):
          if len(scores[a]) > 1:
            score_lim = max(scores[a]) - close_second_best # Allow one mismatch or two consec
            idx = [i for i, score in enumerate(scores[a]) if score > score_lim]
            contacts[a] = [contacts[a][i] for i in idx]
            scores[a] = [scores[a][i] for i in idx]
            
        for a, b in ((0,1), (2,3)):
          if (len(scores[a]) > 1) or (len(scores[b]) > 1):
            n_best_a = scores[a].count(max_score) # Zero is best/perfect score in end-to-end mode, 2 * seq len in local
            n_best_b = scores[b].count(max_score)
            
            if n_best_a * n_best_b == 1: # Only one perfect pair
              i = scores[a].index(max_score)
              j = scores[b].index(max_score)
              
              contacts[a] = [contacts[a][i]]
              contacts[b] = [contacts[b][j]]
            
            else:
              for end in (a, b):
                if len(scores[end]) > 1:
                  chr_name1, start1 = contacts[end][0][0][:2]
                  chr_name2, start2 = contacts[end][1][0][:2]
 
                  if (chr_name1 == chr_name2) and (abs(start2-start1) < CLOSE_AMBIG): # For position ambiguous v. close in cis almost certainly correct
                    i = scores[end].index(max(scores[end]))
                    contacts[end] = [contacts[end][i]]
              
      n_pairs += 1
      iid = int(_id)     
      
      n0 = len(contacts[0])
      n1 = len(contacts[1])
      n2 = len(contacts[2])
      n3 = len(contacts[3])
      
      if (n0+n2) * (n1+n3) == 0: # One end or both ends have no mappings
        n_unmapped += 1
        continue
      
      elif min(n0, n1, n2, n3) == 0: # Not all ends accounted for 
        n_hybrid_end_missing += 1
        continue      

      pairs = []
      for end_a, end_b in ((0,1), (2,3), (0,3), (2,1)): # A:A, B:B, A:B, B:A genome pairings
        for i, (ncc_a, score_a) in enumerate(contacts[end_a]):
          for j, (ncc_b, score_b) in enumerate(contacts[end_b]):
            pairs.append((score_a + score_b, i, j, ncc_a, ncc_b))            
      
      pairs.sort(reverse=True)
      best_score = pairs[0][0]
      score_tol = BOWTIE_MAX_AMBIG_SCORE_TOL
      
      for score, i, j, ncc_a, ncc_b in pairs:
        if score == best_score and (best_score >= really_bad_score):
          chr1 = ncc_a[0]
          chr2 = ncc_b[0]
          
          if ('.' in chr1) and ('.' in chr2):
            chr1, gen1 = chr1.split('.')
            chr2, gen2 = chr2.split('.')
          
            if chr1 == chr2 and (gen1 != gen2):
              score_tol = 3 * BOWTIE_MAX_AMBIG_SCORE_TOL # Stricter for homologous chromosomes
              break
      
      if best_score < really_bad_score: # Nothing any good
        n_hybrid_poor += 1
        continue

      else:
        pairs = [x for x in pairs if x[0] >= best_score-BOWTIE_MAX_AMBIG_SCORE_TOL]
                
      ambig_code = float(len(pairs))
      is_pos_ambig = False
      
      for score, i, j, ncc_a, ncc_b in pairs:
        _write_ncc_line(ncc_a, ncc_b, ambig_code, iid, write_pair)
        ambig_code = 0.0
        if i > 0 or j > 0:
          is_pos_ambig = True
        
      if is_pos_ambig:
        n_ambig += 1
          
      elif len(pairs) > 1: # Ambiguity only due to hybrid genome
        n_hybrid_ambig += 1
      
      else:
        n_unambig += 1

         
  n_primary_strand = [int(x) for x in n_primary_strand]
  n_strand = [int(x) for x in n_strand]
  ncc_file_obj.close()

  stats = [('end_1A_alignments', n_map[0]),
           ('end_2A_alignments', n_map[1]),
           ('end_1B_alignments', n_map[2]),
           ('end_2B_alignments', n_map[3]),
           ('unpaired_ends', (n_unpaired, n_pairs)),
           ('hybrid_all_poor', (n_hybrid_poor, n_pairs)),
           ('hybrid_end_missing', (n_hybrid_end_missing, n_pairs)),
           ('unmapped_end', (n_unmapped, n_pairs)),
           ('unique', (n_unambig, n_pairs)),
           ('genome_ambiguous', (n_hybrid_ambig, n_pairs)),
           ('position_ambiguous', (n_ambig, n_pairs)),
           ('total_contacts', n_pairs)]

  log_report('pair', stats, {'primary_strand':list(zip(n_primary_strand, n_strand))})

  move(paired_ncc_file_name_temp, paired_ncc_file_name)

  return paired_ncc_file_name


def pair_mapped_seqs(sam_file1, sam_file2, chromo_names, file_root,
                     ambig=True, unique_map=False, max_unmap=0):

  paired_ncc_file_name = tag_file_name(file_root, 'pair', '.ncc')
  paired_ncc_file_name_temp = paired_ncc_file_name + TEMP_EXT

  if INTERRUPTED and os.path.exists(paired_ncc_file_name) and not os.path.exists(paired_ncc_file_name_temp):
    return paired_ncc_file_name
   
  ncc_file_obj = open(paired_ncc_file_name_temp, 'w')
  write_pair = ncc_file_obj.write

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
  n_unself = 0
  n_strand1 = 0
  n_strand2 = 0
  n_primary_strand1 = 0
  n_primary_strand2 = 0
  
  max_score = 0
  zero_ord = QUAL_ZERO_ORDS['phred33']
  really_bad_pair_score = 2 * (max_score - 2 * (BOWTIE_MAX_AMBIG_SCORE_TOL+1) )
  close_second_best = 2 * BOWTIE_MAX_AMBIG_SCORE_TOL - 1
  
  # Go through same files and pair based on matching id
  # Write out any to-many mapings to ambiguous

  # Skip to end of headers

  line1 = _skip_sam_header(readline1)
  line2 = _skip_sam_header(readline2)

  # Process data lines

  id1 = line1[:ID_LEN]
  id2 = line2[:ID_LEN]
 
  while id1 and id2:
    if id1 < id2:
      n_unpaired += 1
      line1 = readline1()
      id1 = line1[:ID_LEN]
      continue

    if id2 < id1:
      n_unpaired += 1
      line2 = readline2()
      id2 = line2[:ID_LEN]
      continue

    _id = id1

    contact_a = []
    contact_b = []
    scores_a = []
    scores_b = []

    if id1 != id2:
      raise Exception('Pair mismatch')
    
    j = 0
   
    while id1 == _id:
      
      n_map1 += 1
      data_a = line1.split('\t')
      chr_a = data_a[2]
      
      if chr_a != '*':
        revcomp_a = int(data_a[1]) & 0x10
        start_a = int(data_a[3])
        seq_a = data_a[9]
        qual_a = data_a[10]
        nbp = len(seq_a)
        end_a = start_a + nbp
       
        name_a = chromo_names.get(chr_a, chr_a)
        score_a, var = SCORE_TAG_SEARCH(line1).groups()             
        score_a = int(score_a) # when --local : -nbp*2

        if revcomp_a and var[-1] == '0' and var[-2] in 'GCAT': # Ignore substitutions at the end e.g "C0"; add back subtracted score
          var = var[:-2]
        
          if seq_a[-1] != 'N': # No penalty for N's
            q = min(ord(qual_a[-1]) - zero_ord, 40.0)
            mp1 = 2 + floor(4*q/40.0)  # MX = 6, MN = 2.
            score_a = min(score_a+mp1, max_score)
            end_a -= 1

        elif var[0] == '0' and var[1] in 'GCAT':
          var = var[2:]
          
          if seq_a[0] != 'N':
            q = min(ord(qual_a[0]) - zero_ord, 40.0)
            mp2 = 2 + floor(4*q/40.0)  
            score_a = min(score_a+mp2, max_score)
            start_a += 1

        if revcomp_a:
          start_a, end_a = end_a, start_a  # The sequencing read started from the other end

        ncc_a = (name_a, start_a, end_a, 0, 0, '-' if revcomp_a else '+') # Strand info kept because ends can be different for replicate reads, no Re fragment positions, yet
        contact_a.append((ncc_a, score_a))
        scores_a.append(score_a)
 
        if j == 0:
          if not revcomp_a:
            n_primary_strand1 += 1
 
          n_strand1 += 1

      line1 = readline1()
      id1 = line1[:ID_LEN]
      j += 1
    
    j = 0
    while id2 == _id:
      n_map2 += 1
      data_b = line2.split('\t')
      
      try:
        chr_b = data_b[2]
      except Exception as err:
        warn(line1)
        warn(line2)
        raise err
      
      if chr_b != '*':
        revcomp_b = int(data_b[1]) & 0x10
        start_b = int(data_b[3])
        seq_b = data_b[9]
        qual_b = data_b[10]
        nbp = len(seq_b)
        end_b = start_b + nbp

        name_b = chromo_names.get(chr_b, chr_b)
        score_b, var = SCORE_TAG_SEARCH(line2).groups()             
        score_b = int(score_b) # when --local : -nbp*2

        if revcomp_b and var[-1] == '0' and var[-2] in 'GCAT': # Ignore substitutions at the end e.g "C0"; add back subtracted score
          var = var[:-2]
        
          if seq_b[-1] != 'N': # No penalty for N's
            q = min(ord(qual_b[-1]) - zero_ord, 40.0)
            mp1 = 2 + floor(4*q/40.0)  # MX = 6, MN = 2.
            score_b = min(score_b+mp1, max_score)
            end_b -= 1

        elif var[0] == '0' and var[1] in 'GCAT':
          var = var[2:]
          
          if seq_b[0] != 'N':
            q = min(ord(qual_b[0]) - zero_ord, 40.0)
            mp2 = 2 + floor(4*q/40.0)  
            score_b = min(score_b+mp2, max_score)
            start_b += 1

        if revcomp_b:
          start_b, end_b = end_b, start_b  # The sequencing read started from the other end

        ncc_b = (name_b, start_b, end_b, 0, 0, '-' if revcomp_b else '+')
        contact_b.append((ncc_b, score_b))
        scores_b.append(score_b)
 
        if j == 0:
          if not revcomp_b:
            n_primary_strand2 += 1
 
          n_strand2 += 1

      line2 = readline2()
      id2 = line2[:ID_LEN]
      j += 1

    if not unique_map: # Resolve some perfect vs non-perfect positional ambiguity
      # Assess ambiguous ends separately 
      na = len(scores_a)
      nb = len(scores_b)
      
      # Consider only fairly close second best mappings
      if na > 1:
        score_lim = max(scores_a) - close_second_best # Allow one mismatch or two consec
        idx = [i for i, score in enumerate(scores_a) if score > score_lim] # Mismatches negative; generally best score is 0
        contact_a = [contact_a[i] for i in idx]
        scores_a = [scores_a[i] for i in idx]
        na = len(scores_a)

      if nb > 1:
        score_lim = max(scores_b) - close_second_best # Allow one mismatch or two consec
        idx = [i for i, score in enumerate(scores_b) if score > score_lim]
        contact_b = [contact_b[i] for i in idx]
        scores_b = [scores_b[i] for i in idx]
        nb = len(scores_b)
            
      if 0: # na and nb and (max(na, nb) > 1):
        n_best_a = scores_a.count(max_score) # Zero is best/perfect score in end-to-end mode but not local mode
        n_best_b = scores_b.count(max_score)

        if n_best_a * n_best_b == 1: # Only one perfect pair
          i = scores_a.index(max_score)
          j = scores_b.index(max_score)

          contact_a = [contact_a[i]]
          contact_b = [contact_b[j]]
      
      else:
        # Second best of little consequence if very close in cis
        if len(scores_a) == 2:
          chr_name1, start1 = contact_a[0][0][:2]
          chr_name2, start2 = contact_a[1][0][:2]
 
          if (chr_name1 == chr_name2) and (abs(start2-start1) < CLOSE_AMBIG):
            i = scores_a.index(max(scores_a))
            contact_a = [contact_a[i]]

        if len(scores_b) == 2:
          chr_name1, start1 = contact_b[0][0][:2]
          chr_name2, start2 = contact_b[1][0][:2]
 
          if (chr_name1 == chr_name2) and (abs(start2-start1) < CLOSE_AMBIG):
            i = scores_b.index(max(scores_b))
            contact_b = [contact_b[i]]        
    
    iid = int(_id)
    n_pairs += 1 # Read pairs
    n_a = len(contact_a)
    n_b = len(contact_b)
    n_poss = n_a * n_b # And mapping pairs (maybe ambig)
    
    if n_poss == 0: # Either end not mapped
      n_unmapped += 1
      continue

    elif n_poss == 1: # No multi-position mapping
      n_unambig += 1
      ncc_a, score_a = contact_a[0]
      ncc_b, score_b = contact_b[0]
      _write_ncc_line(ncc_a, ncc_b, 1.0, iid, write_pair)
    
    elif unique_map: # No resoultion of ambiguous ends
      n_ambig += 1
      
      if ambig:
        ambig_code = float(n_poss)
        for ncc_a, score_a in contact_a:
          for ncc_b, score_b in contact_b:
            _write_ncc_line(ncc_a, ncc_b, ambig_code, iid, write_pair) # Write all pair combinations where ambigous
            ambig_code = 0.0   
      
    else:
      # Assess ambiguous paired ends 
     
      pairs = []
      for ncc_a, score_a in contact_a:
        for ncc_b, score_b in contact_b:
          pairs.append((score_a + score_b, ncc_a, ncc_b))
          
      pairs.sort(reverse=True)
      best_score = pairs[0][0]
      
      if best_score < really_bad_pair_score: # Nothing any good; generally a genome build issue
        n_unself += 1
        continue
        
      pairs = [x for x in pairs if x[0] >= best_score-BOWTIE_MAX_AMBIG_SCORE_TOL] 
      n_poss = len(pairs)
      
      if n_poss > 1:
        n_ambig += 1
        
        if ambig:
          ambig_code = float(n_poss)
 
          for score_ab, ncc_a, ncc_b in pairs:
            _write_ncc_line(ncc_a, ncc_b, ambig_code, iid, write_pair) # Write all pair combinations where ambigous
            ambig_code = 0.0
        
      else:
        n_unambig += 1
        score_ab, ncc_a, ncc_b = pairs[0]
        _write_ncc_line(ncc_a, ncc_b, 1.0, iid, write_pair)

  ncc_file_obj.close()

  stats = [('end_1_alignments', n_map1),
           ('end_2_alignments', n_map2),
           ('ambig_all_poor', n_unself),
           ('unpaired_ends', n_unpaired),
           ('unmapped_end', (n_unmapped, n_pairs)),
           ('unique', (n_unambig, n_pairs)),
           ('ambiguous', (n_ambig, n_pairs)),
           ('total_contacts', n_pairs)]

  log_report('pair', stats, {'primary_strand':((n_primary_strand1, n_strand1),
                                               (n_primary_strand2, n_strand2))})

  move(paired_ncc_file_name_temp, paired_ncc_file_name)

  return paired_ncc_file_name


def map_reads(fastq_file, genome_index, align_exe, num_cpu, ambig, qual_scheme, job):

  sam_file_path = tag_file_name(fastq_file, 'map%d' % job, '.sam')
  sam_file_path_temp = sam_file_path + TEMP_EXT

  if INTERRUPTED and os.path.exists(sam_file_path) and not os.path.exists(sam_file_path_temp):
    return sam_file_path

  if os.path.exists(sam_file_path) and not os.path.exists(sam_file_path_temp):
    return sam_file_path

  patt_1 = re.compile('(\d+) reads; of these:')
  patt_2 = re.compile('(\d+) \(.+\) aligned exactly 1 time')
  patt_3 = re.compile('(\d+) \(.+\) aligned 0 times')
  patt_4 = re.compile('(\d+) \(.+\) aligned >1 times')
  
  ##### Micro-C keep unmapped reads ends for another go...
  
  cmd_args = [align_exe,
              '-D', '20', '-R', '3', '-N', '0',  '-L', '20',  '-i', 'S,1,0.5', # similar to very-sensitive
              '-x', genome_index,
              '-k', '2',
              '--reorder',
              '--score-min', 'L,-0.6,-0.6',
              '-p', str(num_cpu),
              '-U', fastq_file,
              '--np', '0', # Penalty for N's
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


def clip_reads(fastq_file, file_root, junct_seq, replaced_seq, qual_scheme, min_qual,
               max_reads_in=None, adapt_seqs=None, trim_5=0, trim_3=0,
               is_second=False, min_len=MIN_READ_LEN):
  """
  Clips reads at ligation junctions, removes anbigous ends and discards very short reads
  """
  
  job = 2 if is_second else 1
  tag = 'reads2_clipped' if is_second else 'reads1_clipped'
  adapt_seqs = adapt_seqs or []
  
  #if is_second and adapt_seqs: # Reverse complement
  #  from string import maketrans
  #  trans_table = maketrans('ATGC', 'TACG')
  #  adapt_seqs += [seq.translate(trans_table)[::-1] for seq in adapt_seqs]
        
  sam_file_path = tag_file_name(file_root, '%s_map%d' % (tag, job), '.sam')
  sam_file_path_temp = sam_file_path + TEMP_EXT
  
  clipped_file = tag_file_name(file_root, tag, '.fastq')
  clipped_file_temp = clipped_file + TEMP_EXT
  
  if INTERRUPTED and os.path.exists(sam_file_path) and not os.path.exists(sam_file_path_temp):
    return clipped_file

  if INTERRUPTED and os.path.exists(sam_file_path) and not os.path.exists(sam_file_path_temp):
    return clipped_file

  if INTERRUPTED and os.path.exists(clipped_file) and not os.path.exists(clipped_file_temp):
    return clipped_file
  
  in_file_obj = open_file_r(fastq_file, complete=not bool(max_reads_in))
  
  if junct_seq:
    n_rep = len(replaced_seq)
    n_junc = len(junct_seq)
  
  else:
    n_rep = 0
    n_junc = 0
  
  re1_pos_counts = defaultdict(int)
  junct_pos_counts = defaultdict(int)
  n_reads = 0
  n_jclip = 0
  n_qclip = 0
  n_short = 0
  n_adapt = 0
  mean_len = 0
  zero_ord = QUAL_ZERO_ORDS[qual_scheme]

  out_file_obj = open(clipped_file_temp, 'w', READ_BUFFER)
  write = out_file_obj.write
  readline = in_file_obj.readline

  line1 = readline()
  while line1[0] != '@':
    line1 = readline()
  
  end = -1 - trim_3
  
  adapt_list = [(adapt_seq, adapt_seq[:MIN_ADAPT_OVERLAP], len(adapt_seq)) for adapt_seq in adapt_seqs]
  
  while line1:
    n_reads += 1
    line2 = readline()[trim_5:end]
    line3 = readline()
    line4 = readline()[trim_5:end]
    
    for adapt_seq, min_adapt, alen in adapt_list:
    
      if min_adapt in line2:
        i = line2.index(min_adapt)
        adapt_end = line2[i:i+alen]
        
        if adapt_end in adapt_seq:
          line2 = line2[:i]
          line4 = line4[:i]
          n_adapt += 1
        
        else:
          n_mismatch = 0
          for j, bp in enumerate(adapt_end):
            if bp == 'N':
              continue
            
            elif bp != adapt_seq[j]:
              n_mismatch += 1
          
          if n_mismatch < 2:
            line2 = line2[:i]
            line4 = line4[:i]
            n_adapt += 1
          
    if junct_seq:
      if junct_seq in line2:
        n_jclip += 1
        i = line2.index(junct_seq)
        junct_pos_counts[i] += 1
        line2 = line2[:i]
 
        if replaced_seq in line2:
          pos = line2.index(replaced_seq)
          re1_pos_counts[pos] += 1
 
        line2 = line2 + replaced_seq
        line4 = line4[:i+n_rep]
 
      elif replaced_seq in line2:
        pos = line2.index(replaced_seq)
        re1_pos_counts[pos] += 1
      
    q = 0
    while line2 and line2[-1] == 'N':
      q = 1
      line2 = line2[:-1]
      line4 = line4[:-1]

    while line4 and (ord(line4[-1]) - zero_ord) < min_qual:
      q = 1
      line2 = line2[:-1]
      line4 = line4[:-1]

    while line4 and (ord(line4[0]) - zero_ord) < min_qual:
      q = 1
      line2 = line2[1:]
      line4 = line4[1:]
    
    for i, qs in enumerate(line4):
      if (ord(qs) - zero_ord) < min_qual:
        line2 = line2[:i] + 'N' + line2[i+1:]
    
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
 
    if max_reads_in and n_reads >= max_reads_in:
      break
    
    line1 = readline()

  if n_reads:
    mean_len /= float(n_reads)
    pos_max = max([max(re1_pos_counts.keys() or [len(line2)]), max(junct_pos_counts.keys() or [len(line2)])])
    re1_pos = list(range(pos_max))
    re1_pos_hist = [100.0 * re1_pos_counts[i]/float(n_reads) for i in re1_pos] 
    junc_pos_hist = [100.0 * junct_pos_counts[i]/float(n_reads) for i in re1_pos] 
    
  else:
    re1_pos = [0]
    re1_pos_hist = [0]
    junc_pos_hist = [0]
  
  in_file_obj.close()
  
  stats = [('input_reads',n_reads),
           ('junction_clipped',(n_jclip, n_reads)),
           ('quaility_clipped',(n_qclip, n_reads)),
           ('too_short',(n_short, n_reads)),
           ('mean_length',mean_len)]
  
  if adapt_seqs:
    stats.insert(3, ("3'_adapter",(n_adapt, n_reads)) )
  
  if is_second:
    stat_key = 'clip_2'
    re1_key = 're1_pos_2'
  else:
    stat_key = 'clip_1'
    re1_key = 're1_pos_1'  
    
  log_report(stat_key, stats, {re1_key:(re1_pos_hist, junc_pos_hist, re1_pos)})

  move(clipped_file_temp, clipped_file)

  return clipped_file

 
def _fragment_fasta(seq, contig, fasta_file_name, frag_size, frag_off):

  info('  .. fragmenting {:,} bp sequence for contig {}'.format(len(seq), contig))

  if os.path.exists(fasta_file_name):
    return 0
   
  fasta_fmt = '>%s__%d__%d__%d__%d\n%s\n'

  n = len(seq)
  out_file_obj = open(fasta_file_name, 'w')
  write = out_file_obj.write

  i = 0
  k = 0
  while i < n:

    j = i+frag_size
    sub_seq = seq[i:j]
    a = b = 0

    while sub_seq and sub_seq[0] == 'N':
      a += 1
      sub_seq = sub_seq[1:]

    while sub_seq and sub_seq[-1] == 'N':
      b += 1
      sub_seq = sub_seq[:-1]

    if len(sub_seq) > MIN_READ_LEN:
      fasta_line = fasta_fmt % (contig, i, j, a, b, sub_seq)
      write(fasta_line)
      k += 1

    i += frag_off

  out_file_obj.close()
  
  return k
  
  
def get_chromo_re_fragments(sequence, re_site, cut_offset):

  frag_data = []
  frag_data_append = frag_data.append
  re_site = re_site.replace('^', '')
  re_site = re_site.replace('_', '')

  if 'N' in re_site:
    re_site = re_site.replace('N', '[AGCT]') # Could allow more ambiguity codes, if feeling generous
    frags = re.split(re_site, sequence)
  else:
    frags = sequence.split(re_site)

  site_len = len(re_site)
  offset_start = site_len - cut_offset
  frag_start = offset_start

  for frag_seq in frags:
    n = len(frag_seq)

    start = frag_start - offset_start
    end = frag_start + n + cut_offset # Offset is cut site, which can be in subseq removed during split

    if not n:
      frag_start += site_len
      continue

    frag_data_append([start, end-1, 1]) # Mappability set later
    frag_start += n + site_len

  frag_data[-1][1] = len(sequence)-1 # Last position is seq end, not cut site

  return frag_data


def get_re_seq(re_name):

  site = RE_SITES[re_name]
  return site.replace('^', '').replace('_','')


def get_chromo_star_pos(sequence, re_sites):
  
  seq_pos = []
  append = seq_pos.append
  
  for re_site in re_sites:
    cut_offset = re_site.replace('_', '').index('^') # From start
    re_site = re_site.replace('^', '').replace('_', '')
    site_len = len(re_site)
    cut_offset = site_len - cut_offset # From end

    if 'N' in re_site:
      re_site = re_site.replace('N', '[AGCT]')
      frags = re.split(re_site, sequence)
    
    else:
      frags = sequence.split(re_site)
    
    frag_start_pos = cut_offset
    for frag_seq in frags:
      n = len(frag_seq)
       
      if not n: # Sites immediately adjascent
        frag_start_pos += site_len
        continue
        
      cut_pos = frag_start_pos - cut_offset # Start of next fragment is after cut site
      append(cut_pos)
        
      frag_start_pos += n + site_len
    
    seq_pos.pop() # Last pos was not an internal site

  return np.array(sorted(seq_pos), int)
  
  
def get_genome_star_sites(re_name, genome_fastas1, genome_fastas2, chromo_name_dict, hom_chromo_dict):
  
  star_sites = STAR_SITES[re_name]
  star_dict = {}
  msg = 'Calculating genome star sites for %s'
  info(msg % re_name)
      
  msg = 'Calculating %s star locations for %s contig %s'
  
  for is_second, fasta_files in enumerate((genome_fastas1, genome_fastas2)):
    if not fasta_files:
      continue
       
    for fasta_file in fasta_files:
      file_obj = open_file_r(fasta_file)
      contig = None
      seq = ''

      for fasta_line in file_obj:
        if fasta_line[0] == '>':
          if contig and seq: # add previous
            chromo = chromo_name_dict.get(contig, contig)
            
            if is_second:
              hybrid_name = chromo + '.b'
            else:
              hybrid_name = chromo + '.a'
 
            if hybrid_name in hom_chromo_dict:
              chr_name = hybrid_name
            else:
              chr_name = chromo
            
            info(msg % (re_name, chr_name, contig))
            star_dict[chr_name] = get_chromo_star_pos(seq.upper(), star_sites)

          contig = fasta_line[1:].split()[0]
          seq = ''

        else:
          seq += fasta_line[:-1]

      file_obj.close()

      if contig and seq:
        chromo = chromo_name_dict.get(contig, contig)
            
        if is_second:
          hybrid_name = chromo + '.b'
        else:
          hybrid_name = chromo + '.a'
 
        if hybrid_name in hom_chromo_dict:
          chr_name = hybrid_name
        else:
          chr_name = chromo
              
        info(msg % (re_name, chr_name, contig))
        star_dict[chr_name] = get_chromo_star_pos(seq.upper(), star_sites)
      
  return star_dict

  
def check_re_frag_file(genome_index, re_name, genome_fastas, chromo_name_dict,
                       align_exe, num_cpu, remap=False):

  re_site_file = get_re_file_name(genome_index, re_name)

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
  
  frag_data = {}
  sources = {}

  msg = 'Calculating %s fragment locations for %s contig %s'
  n_bases = 0

  for fasta_file in genome_fastas:
    file_obj = open_file_r(fasta_file)
    contig = None
    seq = ''

    for fasta_line in file_obj:
      if fasta_line[0] == '>':
        if contig and seq: # add previous
          chromo = chromo_name_dict.get(contig, contig)
          info(msg % (re_name, chromo, contig))
          n_bases += len(seq)
          frag_data[contig] = get_chromo_re_fragments(seq.upper(), re_site, cut_pos)

        contig = fasta_line[1:].split()[0]
        seq = ''

      else:
        seq += fasta_line[:-1]

    file_obj.close()

    if contig and seq:
      chromo = chromo_name_dict.get(contig, contig)
      info(msg % (re_name, chromo, contig))
      n_bases += len(seq)
      frag_data[contig] = get_chromo_re_fragments(seq.upper(), re_site, cut_pos)

  # Write out RE file
  
  contigs = sorted(frag_data)
  
  out_file_obj = open(re_site_file, 'w')
  head = 'Source\tContig\tFrag_Start\tFrag_End\tMappability\n'
  out_file_obj.write(head)

  for contig in contigs:
    source = chromo_name_dict.get(contig, contig)

    for start, end, mappability in frag_data[contig]:
      frag_line = '%s\t%s\t%d\t%d\t%d\n' % (source, contig, start, end, mappability)
      out_file_obj.write(frag_line)

  out_file_obj.close()

  return re_site_file


def read_re_frag_file(file_path, is_second, homo_chromo_dict):

  frag_dict = {}
  chromo_contigs = {}
  npy_cache = file_path + '.npz'

  if os.path.exists(npy_cache):
    in_data_dict = np.load(npy_cache)

    for key in in_data_dict:
      name, contig = key.split()
      frag_dict[name] = in_data_dict[key]

  else:
    file_obj = open_file_r(file_path)
    file_obj.readline()

    prev_contig = None
    for line in file_obj:
      name, contig, start, end, mappability = line.split()

      if contig != prev_contig:
        frag_dict[name] = []
        append = frag_dict[name].append
        chromo_contigs[name] = contig

      append((int(start), int(end), int(mappability)))
      prev_contig = contig

    file_obj.close()

    kw_args = {}
    for name, contig in chromo_contigs.items():
      key = '%s %s' % (name, contig)
      kw_args[key] = np.array(frag_dict[name])

    try:
      np.savez(npy_cache, **kw_args)
      
    except Exception:
      pass # Not to worry if cache fails: it's only a speedup

  end_dict = {}
  out_frag_dict = {} # Carries hybrid chromosome names, as appropriate
  
  for name in frag_dict:
    if is_second:
      hybrid_name = name + '.b'
    else:
      hybrid_name = name + '.a'
    
    if hybrid_name in homo_chromo_dict:
      chr_name = hybrid_name
    else:
      chr_name = name
    
    frag_regions = np.array(frag_dict[name])
    out_frag_dict[chr_name] = frag_regions
    end_dict[chr_name] = frag_regions[:,1]

  return out_frag_dict, end_dict


def uncompress_gz_file(file_name):

  if file_name.endswith('.gz'):
    in_file_obj = gzip.open(file_name, 'rb')

    file_name = file_name[:-3]
    out_file_obj = open(file_name, 'w')
    write = out_file_obj.write
    
    for line in in_file_obj:
      write(line)

    # Faster alternative, given sufficient memory:
    # out_file_obj.write(infile_obj.read())

    in_file_obj.close()
    out_file_obj.close()

  return file_name


def index_genome(base_name, file_names, output_dir, indexer_exe='bowtie2-build',
                 table_size=10, quiet=True, pack=True, num_cpu=2):

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

  cmd_args += ['-t', str(table_size), '--threads', str(num_cpu), fasta_file_str, base_name]
  #cmd_args += ['-t', str(table_size), fasta_file_str, base_name]

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


def check_re_site(seq, critical=True):
   
  is_ok = True
  msgs = []
   
  if set(seq) - set('^_GCATN'):
    msg = 'Restriction site specification (%s) may only contain characters "G", "C", "A", "T", "N", "_", or "^".' % seq
    is_ok = False
    msgs.append(msg)

  if '^' not in seq:
    msg = 'Restriction site specification (%s) does not contain a "^" cut point.' % seq
    is_ok = False
    msgs.append(msg)

  if seq.count('^') > 1:
    if '_' in seq:
      msg = 'Restriction site specification (%s) may only contain one "^" 5\' cut point.' % seq

    else:
      msg = 'Restriction site specification (%s) may only contain one "^" 5\' cut point, use "_" for 3\'.' % seq

    is_ok = False
    msgs.append(msg)

  if seq.count('_') > 1:
    msg = 'Restriction site specification (%s) may only contain one "_" 3\' cut point.' % seq
    is_ok = False
    msgs.append(msg)
  
  msg = ' '.join(msgs)
  
  if critical and not is_ok:
    fatal(msg)

  return is_ok, msg


def check_fastq_file(file_path, critical=True):

  is_ok = True
  msg = ''
  file_obj = open_file_r(file_path, complete=False)
  
  lines = file_obj.readlines(FASTQ_READ_CHUNK)
  lines = [l for l in lines if l.rstrip()]

  for i in range(0, len(lines), 4):
    if not (lines[i][0] == '@') and (lines[i+2][0] == '+'):
      msg = 'File "%s" does not appear to be in FASTQ format' % file_path
      is_ok = False
  
  file_obj.close()
  
  if critical and not is_ok:
    fatal(msg)

  return is_ok, msg


def is_genome_indexed(file_path, file_exts=('.1.bt2','.1.bt2l')):

  for file_ext in file_exts:
    if os.path.exists(file_path + file_ext):
      return True

  return False


def check_index_file(file_path, sub_files=('.1', '.2', '.3', '.4', '.rev.1', '.rev.2'), critical=True):
     
  msg = ''
  is_ok = True

  if os.path.exists(file_path + '.1.bt2l'):
    file_ext = '.bt2l' # All build files should be long
  else:
    file_ext = '.bt2'

  for sub_file in sub_files:
    full_path = file_path + sub_file + file_ext
    is_ok, msg = check_regular_file(full_path)
    
    if not is_ok:
      msg = 'Genome index error. ' + msg    
      break

  if critical and not is_ok:
    fatal(msg)

  return is_ok, msg


def check_file_extension(file_path, file_ext):

  file_ext = file_ext.lower()

  if not file_ext.startswith('.'):
    file_ext = '.' + file_ext

  if not file_path.lower().endswith(file_ext):
    file_path = os.path.splitext(file_path)[0] + file_ext

  return file_path


def check_regular_file(file_path, critical=False):

  msg = ''
  is_ok = True

  if not os.path.exists(file_path):
    msg = 'File "%s" does not exist' % file_path
    is_ok = False

  elif not os.path.isfile(file_path):
    msg = 'Location "%s" is not a regular file' % file_path
    is_ok = False

  elif os.stat(file_path).st_size == 0:
    msg = 'File "%s" is of zero size ' % file_path
    is_ok = False

  elif not os.access(file_path, os.R_OK):
    msg = 'File "%s" is not readable' % file_path
    is_ok = False

  if critical and not is_ok:
    fatal(msg)

  return is_ok, msg


def _write_log_lines(lines, verbose=VERBOSE):

  if LOG_FILE_PATH:
    with open(LOG_FILE_PATH, 'a') as file_obj:
      file_obj.write('\n'.join(lines) + '\n')

  if verbose:
    for line in lines:
      print(line)


def info(msg, prefix='INFO'):

  line = '%s: %s' % (prefix, msg)
  _write_log_lines([line])


def warn(msg, prefix='WARNING'):

  line = '%s: %s' % (prefix, msg)
  _write_log_lines([line])


def fatal(msg, prefix='%s FAILURE' % PROG_NAME):

  lines = ['%s: %s' % (prefix, msg),
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

  file_obj = open_file_r(file_path, complete=False)

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


def get_re_file_name(genome_index, re_name):

  dir_name, index_file = os.path.split(genome_index)
  re_site_file = 'RE_frag_%s_%s.txt' % (re_name, index_file)
  
  return os.path.join(dir_name, re_site_file) 


def read_chromo_names(chr_name_paths, genome_indices, re_name):
  
  hom_chromo_dict = {}
  name_dicts = [{}, {}]
  contig_dicts = [{}, {}]
  
  for i, genome_index in enumerate(genome_indices):
    
    if chr_name_paths[i]: # Chromosome naming file takes precedence
      file_path = chr_name_paths[i]
      msg = 'Using chromosome names from file %s'
      info(msg % file_path)
      check_regular_file(file_path, critical=True)
      
      with open(file_path, 'rU') as file_obj:
        for line in file_obj:
          line = line.strip()
 
          if not line:
            continue
 
          if line[0] == '#':
            continue
 
          data = line.split()
 
          if len(data) < 2:
            continue
 
          else:
            contig, name = data[:2]
 
            if contig in name_dicts[i]:
              msg ='Chromosome naming file "%s" contains a repeated sequence/contig name: %s'
              fatal(msg % (file_path, contig))
 
            if name in contig_dicts[i]:
              msg ='Chromosome naming file "%s" contains a repeated chromosome/segment name: %s'
              fatal(msg % (file_path, name))
 
            name_dicts[i][contig] = name
            contig_dicts[i][name] = contig
 
      if not name_dicts[i]:
        fatal('Chromosome naming file "%s" contained no usable data: requires whitespace-separated pairs' % file_path)      
    
    elif re_name: # Use any existing RE file
      file_path = get_re_file_name(genome_index, re_name)
      
      if os.path.exists(file_path): # Load names from RE file
        check_regular_file(file_path, critical=True) 
        msg = 'Using chromosome names from RE mapping file %s'
        info(msg % file_path)
        npy_cache = file_path + '.npz'

        if os.path.exists(npy_cache):
          in_data_dict = np.load(npy_cache)

          for key in in_data_dict:
            name, contig = key.split()
            name_dicts[i][contig] = name
            contig_dicts[i][name] = contig

        else:
          with open_file_r(file_path) as file_obj:
            null = file_obj.readline()
            prev_contig = None
            
            for line in file_obj:
              name, contig, start, end, mappability = line.split()

              if contig != prev_contig:
                name_dicts[i][contig] = name
                contig_dicts[i][name] = contig
               
              prev_contig = contig
 
      else:
        msg = 'A chromosome naming file must be specified for genome index (%s) as ' \
              'this required to create its RE mapping file'
        fatal(msg % genome_index)
    else:
      msg = 'A chromosome naming file must be specified for genome index (%s)'
      fatal(msg % genome_index)
  
  orig_name_dicts = [dict(x) for x in name_dicts] 
  
  if contig_dicts[0] and contig_dicts[1]:
    for name in list(contig_dicts[0].keys()):
      if name in contig_dicts[1]:
        contig1 = contig_dicts[0][name]
        contig2 = contig_dicts[1][name]
    
        name_a = name + '.a'
        name_b = name + '.b'
        
        hom_chromo_dict[contig1] = contig2
        hom_chromo_dict[contig2] = contig1

        hom_chromo_dict[name_a] = name_b
        hom_chromo_dict[name_b] = name_a

        name_dicts[0][contig1] = name_a
        name_dicts[1][contig2] = name_b
    
  chr_name_dict = name_dicts[0]
  chr_name_dict.update(name_dicts[1])
  
  base_chromo_name_dict = orig_name_dicts[0]
  base_chromo_name_dict.update( orig_name_dicts[1])
  
  return chr_name_dict, base_chromo_name_dict, hom_chromo_dict
  

def nuc_process(fastq_paths, genome_index, genome_index2, re1, re2=None, chr_names1=None, chr_names2=None,
                sizes=(300,800), min_rep=2, num_cpu=1, num_copies=1, ambig=True, unique_map=False, out_file=None,
                report_file=None, align_exe=None, qual_scheme=None, min_qual=30, g_fastas=None,
                g_fastas2=None, is_pop_data=False, remap=False, reindex=False, keep_files=True, lig_junc=None,
                zip_files=True, sam_format=True, verbose=True, max_reads_in=None, adapt_seqs=None,
                trim_5=0, trim_3=0):
  """
  Main function for command-line operation
  """
  
  from hic_core.nuc_process_report import nuc_process_report
  from hic_core.sc_hic_disambiguate import sc_hic_disambiguate
  from tools.contact_map import contact_map
  from tools.ncc_bin import bin_ncc
  
  genome_indices = [genome_index]
  if genome_index2:
    if not (chr_names1 and chr_names2):
      fatal('Two chromosome naming files must be specified (-c and -c2) for dual genome indices')

    genome_indices.append(genome_index2)
    is_hybrid = True
  else:
    is_hybrid = False 
  
  if (re1 is None) and not chr_names1:
     fatal('A chromosome naming file must be specified (-c) when primary restriction enzyme (re1) is None')
    
  # Files can be empty if just re-indexing etc...
  if not fastq_paths and not (reindex and g_fastas):
    fatal('No FASTQ files specified')

  # Check FASTQ files : two present, exist, right format
  for file_path in fastq_paths:
    check_regular_file(file_path, critical=True)
    check_fastq_file(file_path, critical=True)
    
  if fastq_paths:
    if len(fastq_paths) == 1:
      fatal('Only one FASTQ file specified (exactly two required)')

    if len(fastq_paths) > 2:
      fatal('More than two FASTQ files specified (exactly two required)')

  # Ambiguous data is usually only for single-cell use
  if is_pop_data and ambig:
    warn('Ambiguous option is being used with population Hi-C data')

  if not (0 <= min_qual <= 40):
    fatal('Miniumum FASTQ quality score must be in the range 0-40 (%d specified).' % min_qual)

  # Check genome index files, if not being rebuilt
  for gidx in genome_indices:
    if is_genome_indexed(gidx):
      check_index_file(gidx)
    else:
     dir_name, idx_name = os.path.split(gidx)
   
     if not os.path.exists(dir_name):
       msg = 'Genome index location "{}" does not exist'.format(dir_name)
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
    msg = 'Aligner executable path "%s" not found'
    fatal(msg % align_exe)

  else:
    check_regular_file(align_exe, critical=True)

  # Check reindexing input present
  if reindex and not g_fastas:
    fatal('No genome FASTA files specified for re-indexing')
  
  if re1 in STAR_SITES:
    if not g_fastas:
      fatal('No genome FASTA files specified for checking %s star sites' % re1)
    elif is_hybrid and not g_fastas2:
      fatal('No secondary genome FASTA files specified for checking %s star sites' % re1)
  
  # Must have a preconstructed genome index or FASTA files to build one from
  if not (genome_index or g_fastas):
    fatal('No genome index file of genome FASTA sequences specified')

  # Check genome FASTA files
  if g_fastas:
    for fasta_path in g_fastas:
      check_regular_file(fasta_path, critical=True)

  if g_fastas2:
    for fasta_path in g_fastas2:
      check_regular_file(fasta_path, critical=True)

  # Check restriction enzymes
  if re1 and re1 not in RE_SITES:
    msg = 'Restriction enzyme "%s" not known. Available: %s.' % (re1, ', '.join(sorted(RE_SITES)))
    msg += ' %s can be edited to add further cut-site defintions.' % RE_CONF_FILE
    fatal(msg)
  
  if re1:
    check_re_site(RE_SITES[re1], critical=True)
    re1Seq = get_re_seq(re1)
  
  else:
    re1Seq = None
  
  # TBC: re2Seq is never used
  if re2:
    re2Seq = get_re_seq(re2)

    if re2 not in RE_SITES:
      msg = 'Restriction enzyme "%s" not known. Available: %s.' % (re2, ', '.join(sorted(RE_SITES)))
      msg += ' Enzymes.conf can be edited to add further cut-site defintions.'
      fatal(msg)

    check_re_site(RE_SITES[re2], critical=True)

  else:
    re2Seq = None
    
  # Check and get any chromosome naming files
  chromo_name_dict, base_chromo_name_dict, hom_chromo_dict = read_chromo_names((chr_names1, chr_names2), genome_indices, re1)
  
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

  elif re1:
    lig_junc = get_ligation_junction(RE_SITES[re1])
  
  else:
    lig_junc = None
  
  # Check adapter sequences
  if adapt_seqs:
    adapt_seqs = [x.upper() for x in adapt_seqs]
    for adapt_seq in adapt_seqs:
      if set(adapt_seq) - set('GCAT'):
        msg = 'Adapter sequence may only contain letters G, C, A and T'
        fatal(msg)    
    
  # Get base file name for output
  for file_path in (out_file, report_file):
    if file_path:
      file_root = os.path.splitext(file_path)[0]
      break

  else:
    file_paths = []
    for fastq_path in fastq_paths:
      if fastq_path.lower().endswith('.gz'):
        fastq_path = fastq_path[:-3]

      file_paths.append(fastq_path)

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
    ambig_file = file_root + '_ambig.ncc'

  else:
    ambig_file = None

  if report_file:
    report_file = check_file_extension(report_file, '.pdf')
  else:
    report_file = file_root + '_report.pdf'

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
    
  log_report('command', [('Command used', ' '.join(sys.argv)),('Version', VERSION)])

  general_stats = [
    ('Input FASTQ file 1', os.path.basename(fastq_paths[0])),
    ('Input FASTQ file 2', os.path.basename(fastq_paths[1])),
    ('Output contacts file', os.path.basename(out_file)),
    ('Report file', os.path.basename(report_file)),
    ('Output directory', os.path.dirname(os.path.abspath(out_file))),
    ('Aligner executable', align_exe),
    ('Genome indices', ', '.join(genome_indices)),
    ('FASTQ quality scheme', qual_scheme),
    ('Minimum FASTQ quaility', min_qual),
    ('RE1 site', (re1 or 'None')),
    ('RE2 site', (re2 or 'None')),
    ('Ligation junction', lig_junc),
    ('Adapter sequences', ','.join(adapt_seqs or ['None'])),
    ('Molecule size range', size_str),
    ('Min sequencing repeats', min_rep),
    ('Parallel CPU cores', num_cpu),
    ('End trimming 5\', 3\'', '%d, %d' % (trim_5, trim_3)),
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
    index_genome(base_name, g_fastas, output_dir or '.', indexer_exe, num_cpu=num_cpu) # Latest version of bowtie2 can do parallel index builds (--threads)

  if g_fastas2 and genome_index2 and (reindex or not is_genome_indexed(genome_index2)):
    warn('Indexing secondary genome, this may take some time...')
    output_dir, base_name = os.path.split(genome_index2)
    index_genome(base_name, g_fastas2, output_dir or '.', indexer_exe, num_cpu=num_cpu)
  
  
  # Check for extra RE1 cutting at star sites, for some enzymes (e.g. HindIII)
  
  star_dict = {}
  if re1 in STAR_SITES:
    star_dict = get_genome_star_sites(re1, g_fastas, g_fastas2, base_chromo_name_dict, hom_chromo_dict)
    
  # Create RE fragments file if not present
  
  re1_files = []
  re2_files = []
  
  if re1:
    re1_files = [check_re_frag_file(genome_index, re1, g_fastas, base_chromo_name_dict,
                                    align_exe, num_cpu, remap=remap)]

  if re2:
    re2_files.append(check_re_frag_file(genome_index, re2, g_fastas, base_chromo_name_dict,
                                        align_exe, num_cpu, remap=remap))
  
  if is_hybrid:
    if re1:
      re1_files.append(check_re_frag_file(genome_index2, re1, g_fastas2, base_chromo_name_dict,
                                          align_exe, num_cpu, remap=remap))

    if re2:
      re2_files.append(check_re_frag_file(genome_index2, re2, g_fastas2, base_chromo_name_dict,
                                          align_exe, num_cpu, remap=remap))

  # Clip read seqs at any sequenced ligation junctions  

  info('Clipping FASTQ reads...')
  
  clipped_file1 = clip_reads(fastq_paths[0], intermed_file_root, lig_junc, re1Seq,
                             qual_scheme, min_qual, max_reads_in, adapt_seqs,
                             trim_5, trim_3)
                             
  clipped_file2 = clip_reads(fastq_paths[1], intermed_file_root, lig_junc, re1Seq,
                             qual_scheme, min_qual, max_reads_in, adapt_seqs,
                             trim_5, trim_3, is_second=True)
                             
  # Do the main genome mapping
  info('Mapping FASTQ reads...')
  sam_file1 = map_reads(clipped_file1, genome_index, align_exe, num_cpu, ambig, qual_scheme, 1)
  sam_file2 = map_reads(clipped_file2, genome_index, align_exe, num_cpu, ambig, qual_scheme, 2)

  if is_hybrid:
    info('Mapping FASTQ reads to second genome...')
    sam_file3 = map_reads(clipped_file1, genome_index2, align_exe, num_cpu, ambig, qual_scheme, 3)
    sam_file4 = map_reads(clipped_file2, genome_index2, align_exe, num_cpu, ambig, qual_scheme, 4)
  
  if os.path.exists(clipped_file2) and not keep_files:
    os.unlink(clipped_file1)
    os.unlink(clipped_file2)  
  
  info('Pairing FASTQ reads...')

  if is_hybrid:
    paired_ncc_file = pair_mapped_hybrid_seqs(sam_file1, sam_file2, sam_file3, sam_file4,
                                              chromo_name_dict, intermed_file_root, ambig, unique_map)
  else:
    paired_ncc_file = pair_mapped_seqs(sam_file1, sam_file2, chromo_name_dict,
                                       intermed_file_root, ambig, unique_map)

  # Write SAM, if requested
  if sam_format:
    write_sam_file(paired_ncc_file, sam_file1, sam_file2)

  # Filtering
  info('Filtering mapped sequences...')
   
  filter_output = filter_pairs(paired_ncc_file, re1_files, re2_files, chromo_name_dict,
                               hom_chromo_dict, sizes, keep_files, zip_files, star_dict=star_dict)
  
  filter_ncc_file, fail_file_names = filter_output

  # Write SAM, if requested
  if sam_format:
    write_sam_file(filter_ncc_file, sam_file1, sam_file2)

    for key in fail_file_names:
      write_sam_file(fail_file_names[key], sam_file1, sam_file2)

  # Redundancy reduction and promiscuity checks not relevant for population data
  # where multiple observations from and between the same RE fragments are expected

  if is_pop_data:
    
    if re2:
      nr_ncc_file = None
      
    else:  
      # De-duplication of polulation data only if exact read bp matches and no fragment redundency requirement 
      info('Removing duplicate contacts... ')
      nr_ncc_file = remove_redundancy(filter_ncc_file, keep_files, zip_files, use_re_fragments=False)

    if sam_format and nr_ncc_file:
      write_sam_file(nr_ncc_file, sam_file1, sam_file2)   

    if keep_files:    
      if nr_ncc_file:
        shutil.copyfile(nr_ncc_file, out_file)
      else:
        shutil.copyfile(filter_ncc_file, out_file)
        
      if zip_files and nr_ncc_file:
        nr_ncc_file = compress_file(nr_ncc_file)

    elif nr_ncc_file:
      shutil.move(nr_ncc_file, out_file)
     
    else:
      shutil.move(filter_ncc_file, out_file)
       
  else:
    # Merge duplicates
    info('Removing duplicate contacts... ')
    
    if re1:
      nr_ncc_file = remove_redundancy(filter_ncc_file, keep_files, zip_files, min_rep, is_hybrid=is_hybrid)
    else: # Whene no RE, e.g. with MicroC, redundancy removal is only for exact duplicates, c.f. population data
      nr_ncc_file = remove_redundancy(filter_ncc_file, keep_files, zip_files, use_re_fragments=False, is_hybrid=is_hybrid)
    
    if sam_format:
      write_sam_file(nr_ncc_file, sam_file1, sam_file2)

    # Remove promiscuous ends
    info('Removing pairs with promiscuous ends...')
    clean_ncc_file = remove_promiscuous(nr_ncc_file, num_copies, keep_files, zip_files)

    if sam_format:
      write_sam_file(clean_ncc_file, sam_file1, sam_file2)

    if keep_files:
      shutil.copyfile(clean_ncc_file, out_file)
      
      if zip_files:
        clean_ncc_file = compress_file(clean_ncc_file)

    else:
      os.unlink(nr_ncc_file)
      shutil.move(clean_ncc_file, out_file)
  
  if not keep_files:
    os.unlink(sam_file1)
    os.unlink(sam_file2)
    os.unlink(paired_ncc_file)
    
    if is_hybrid:
      os.unlink(sam_file3)
      os.unlink(sam_file4)
    
    try:
      os.rmdir(intermed_dir)
    
    except OSError as err:     
      warn('Could not remove processing directory {}. Python error: {}'.format(intermed_dir, err))
      
  final_stats = get_ncc_stats(out_file, hom_chromo_dict)
  
  log_report('final', final_stats)
  
  nuc_process_report(STAT_FILE_PATH, report_file)

  if final_stats[0][1] > 1: # Final contacts
    pdf_path = '{}_contact_map.pdf'.format(os.path.splitext(out_file)[0])
    # Contact_map resolution to be determined by chromosome sizes
    
    if is_pop_data:
      npz_path = '{}.npz'.format(os.path.splitext(out_file)[0])
      bin_ncc(out_file, npz_path, bin_size=25.0)
      contact_map([npz_path], pdf_path, bin_size=None, bin_size2=None,
                  no_separate_cis=False, is_single_cell=False)

    else:
      contact_map([out_file], pdf_path, bin_size=None, bin_size2=None,
                  no_separate_cis=False, is_single_cell=True)
  
  if is_hybrid and not is_pop_data:
    sc_hic_disambiguate([out_file])
  
                
  info('Nuc Process all done.')

  return out_file


def main(argv=None):
  from argparse import ArgumentParser

  if argv is None:
    argv = sys.argv[1:]

  epilog = 'Note %s can be edited to add further restriction enzyme cut-site defintions. ' % RE_CONF_FILE
  epilog += 'For further help email tjs23@cam.ac.uk or wb104@cam.ac.uk'

  res = sorted(RE_SITES) + ['None']
  avail_re = 'Available: ' + ', '.join(res)
  avail_quals = 'Available: ' + ', '.join(QUAL_SCHEMES)

  arg_parse = ArgumentParser(prog=PROG_NAME, description=DESCRIPTION,
                             epilog=epilog, prefix_chars='-', add_help=True)

  arg_parse.add_argument(nargs='+', metavar='FASTQ_FILE', dest='i',
                         help='Input paired-read FASTQ files to process. Accepts wildcards' \
                              ' that match paired files. If more than two files are input,' \
                              ' processing will be run in batch mode using the same parameters.')

  arg_parse.add_argument('-g', '--gnome-index', metavar='GENOME_FILE', dest='g',
                         help='Location of genome index files to map sequence reads to without' \
                              ' any file extensions like ".1.b2" etc. A new index will be created' \
                              ' with the name if the index is missing and genome FASTA files are' \
                              ' specified')

  arg_parse.add_argument('-g2', '--gnome-index-2', default=None, metavar='GENOME_FILE_2', dest='g2',
                         help='Location of secondary genome index files for hybrid genomes. A new' \
                              ' index will be created with the name if the index is missing and' \
                              ' genome FASTA files are specified')
  
  arg_parse.add_argument('-cn', '--chromo-names', default=None, metavar='CHROM_NAME_FILE', dest='cn',
                         help='Location of a file containing chromosome names for the genome build:' \
                              ' tab-separated lines mapping sequence/contig names (as appear at the' \
                              ' start of genome FASTA headers) to desired (human readable) chromosome' \
                              ' names. This file is not mandatory if the primary restriction enzyme (-re1)' \
                              ' is specified (i.e. not "None") and a corresponding RE1 mapping file has already' \
                              ' been created for the genome. The naming file may be built automatically from NCBI' \
                              ' genome FASTA files using the supplied "nuc_sequence_names" program')  

  arg_parse.add_argument('-cn2', '--chromo-names-2', default=None, metavar='CHROM_NAME_FILE_2', dest='cn2',
                         help='Location of a file containing chromosome names for a second' \
                              ' hybrid genome build. This file is only mandatory if an RE1 mapping file has not already' \
                              ' been created for the genome. The names in this file must exactly match those' \
                              ' for the other genome build where chromosomes are homologous. The' \
                              ' file may be built automatically from NCBI genome FASTA files using' \
                              ' the supplied "nuc_sequence_names" program')  
  
  arg_parse.add_argument('-re1', default='MboI', choices=res, metavar='ENZYME',
                         help='Primary restriction enzyme (for ligation junctions). ' \
                              'May be set to "None" for MicroC etc., where digestion is not sequence specific. ' \
                              'Default: MboI.' + avail_re)

  arg_parse.add_argument('-re2', choices=res, metavar='ENZYME',
                         help='Secondary restriction enzyme (if used). ' + avail_re)

  arg_parse.add_argument('-s', '--size-range', default='150-1500', metavar='SIZE_RANGE', dest='s',
                         help='Allowed range of sequenced molecule sizes,' \
                              ' e.g. "150-1000", "100,800" or "200" (no maximum)')

  arg_parse.add_argument('-n', '--num-cpu', default=1, metavar='CPU_COUNT', dest='n',
                         type=int, help='Number of CPU cores to use in parallel')

  arg_parse.add_argument('-r', '--redundant-min', default=2, metavar='COUNT', type=int, dest='r',
                         help='Minimum number of redundant sequencing repeats (e.g. PCR amplicons)' \
                              ' required to support a contact')

  arg_parse.add_argument('-o', '--out-contact-file', metavar='NCC_FILE', dest='o',
                         help='Optional output name for NCC format chromosome contact file. This' \
                              ' option will be ignored if more than two paired FASTA files are' \
                              ' input (i.e. for batch mode); automated naming will be used instead.')

  arg_parse.add_argument('-pdf', '--pdf-report-file', metavar='PDF_FILE', dest='pdf',
                         help='Optional output name for PDF format report file. This option will' \
                              ' be ignored if more than two paired FASTA files are input (i.e.' \
                              ' for batch mode); automated naming will be used instead.')

  arg_parse.add_argument('-b', '--bowtie2-path', metavar='EXE_FILE', dest='b',
                         help='Path to bowtie2 (read aligner) executable (will be searched' \
                              ' for if not specified)')

  arg_parse.add_argument('-q', '--qual-scheme', metavar='SCHEME', dest='q',
                         help='Use a specific FASTQ quality scheme (normally not set' \
                              ' and deduced automatically). ' + avail_quals)

  arg_parse.add_argument('-qm', '--qual-min', default=DEFAULT_MIN_QUAL, metavar='MIN_QUALITY', type=int, dest='qm',
                         help='Minimum acceptable FASTQ quality score in range 0-40 for' \
                              ' clipping end of reads. Default: %d' % DEFAULT_MIN_QUAL)

  arg_parse.add_argument('-m', '--map-re-sites', default=False, action='store_true', dest='m',
                         help='Force a re-mapping of genome restriction enzyme sites' \
                              ' (otherwise cached values will be used if present)')

  arg_parse.add_argument('-p', '--population', default=False, action='store_true', dest='p',
                         help='The input data is multi-cell/population Hi-C;' \
                              ' single-cell processing steps are avoided')

  arg_parse.add_argument('-pt', '--pair-end-tags', nargs=2, metavar='PAIRED_READ_TAGS', dest='pt', default=['r_1','r_2'],
                         help='When more than two FASTQ files are input (batch mode),' \
                              ' the subtrings/tags which differ between paired FASTQ' \
                              ' file paths. Default: r_1 r_2')

  arg_parse.add_argument('-x', '--reindex', default=False, action='store_true', dest='x',
                         help='Force a re-indexing of the genome (given appropriate FASTA files)')

  arg_parse.add_argument('-f', '--genome-fastas', nargs='+', metavar='FASTA_FILES', dest='f',
                         help='Specify genome FASTA files for genome index building' \
                              ' (accepts wildcards)')

  arg_parse.add_argument('-f2', '--genome-fastas-2', nargs='+', metavar='FASTA_FILES_2', dest='f2',
                         help='A second set of genome FASTA files for building a second genome' \
                              ' index when using hybrid strain cells (accepts wildcards).')

  arg_parse.add_argument('-a', '--ambiguous', default=False, action='store_true', dest='a',
                         help='Whether to report ambiguously mapped contacts')

  arg_parse.add_argument('-k', '--keep-files', default=False, action='store_true', dest='k',
                         help='Keep any intermediate files (e.g. clipped FASTQ etc).')

  arg_parse.add_argument('-sam', default=False, action='store_true',
                         help='Write paired contacts files to SAM format')

  arg_parse.add_argument('-l', '--ligation-junct', metavar='SEQUENCE', dest='l',
                         help='Seek a specific ligation junction sequence (otherwise this is' \
                              ' guessed from the primary restriction enzyme)')

  arg_parse.add_argument('-z', '--gzip-output', default=False, action='store_true', dest='z',
                         help='GZIP compress any output FASTQ files')

  arg_parse.add_argument('-v', '--verbose', action='store_true', dest='v',
                         help='Display verbose messages to report progress')

  arg_parse.add_argument('-u', '--unique-only', default=False, action='store_true', dest='u',
                         help='Whether to only accept uniquely mapping genome positions and not' \
                              ' attempt to resolve certain classes of ambiguous mapping where a' \
                              ' single perfect match is found.')

  arg_parse.add_argument('-cc', '--chromo-copies', default=0, metavar='GENOME_COPIES', dest='cc',
                         type=int, help='Number of whole-genome copies, e.g. for G2 phase;' \
                              ' Default 1 unless as second genome index is specified' \
                              ' for hybrid samples.')

  arg_parse.add_argument('-lim', '--limit-reads', default=0, metavar='MAX_READS', dest='lim',
                         type=int, help='Limit the number of input reads considered: useful for' \
                         'testing population Hi-C data prior to a lengthy full run')

  arg_parse.add_argument('-5', '--5-prime-trim', default=0, metavar='CLIP_BP', dest='5',
                         type=int, help='Number of base pairs to trim from the 5\' start of all input reads.' \
                         'Default: 0')

  arg_parse.add_argument('-3', '--3-prime-trim', default=0, metavar='CLIP_BP', dest='3',
                         type=int, help='Number of base pairs to trim from the 3\' end of all input reads.' \
                         'Default: 0')

  default_ad_seq = ADAPTER_SEQS[DEFAULT_ADAPTER]
  ad_prests = ', '.join(['%s:%s' % (k, v) for k, v in ADAPTER_SEQS.items()])
  
  arg_parse.add_argument('-ad', '--adapter-seq', nargs='*',  default=default_ad_seq,
                         metavar='ADAPTER_SEQ', dest='ad',
                         help='Adapter sequences to truncate reads at (or blank for none). E.g. %s. ' \
                         'Default: %s (%s)' % (ad_prests, ADAPTER_SEQS[DEFAULT_ADAPTER], DEFAULT_ADAPTER))

  args = vars(arg_parse.parse_args(argv))

  fastq_paths = args['i']
  genome_index = args['g']
  genome_index2 = args['g2']
  re1 = args['re1']
  re2 = args['re2']
  c1 = args['cn']
  c2 = args['cn2']
  sizes = args['s']
  num_cpu = args['n']
  min_rep = args['r']
  ambig = args['a']
  out_file = args['o']
  report_file = args['pdf']
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
  num_copies = args['cc']
  lim_reads = args['lim'] or None
  adapt_seqs = args['ad']
  trim_5 = args['5']
  trim_3 = args['3']
  
  if adapt_seqs and isinstance(adapt_seqs, str):
    adapt_seqs = [adapt_seqs]
  
  if not num_copies:
    if genome_index2:
      num_copies = 4
    else:
      num_copies = 2
  
  if re1.lower() == 'none':
    re1 = None
  
  if sizes:
    sizes = sorted([int(x) for x in re.split('\D+', sizes)])

  if not fastq_paths:
    msg = 'No FASTQ paths specified'
    fatal(msg)

  elif len(fastq_paths) > 2: # Batch mode
    fastq_paths_1, fastq_paths_2 = pair_fastq_files(fastq_paths, pair_tags)
    fastq_pairs = list(zip(fastq_paths_1, fastq_paths_2))

    if out_file or report_file:
      msg = 'Output naming options ignored for batch mode. Automated naming will be used for each FASTQ pair'
      warn(msg)
      out_file = None
      ambig_file = None
      report_file = None
  
  else:
    fastq_pairs = [fastq_paths]
  
  for fastq_path_pair in fastq_pairs:
    nuc_process(fastq_path_pair, genome_index, genome_index2, re1, re2, c1, c2, sizes, min_rep,
                num_cpu, num_copies, ambig, unique_map, out_file, report_file, align_exe,
                qual_scheme, min_qual, g_fastas, g_fastas2, is_pop_data, remap, reindex, keep_files,
                lig_junc, zip_files, sam_format, verbose, lim_reads, adapt_seqs, trim_5, trim_3)

  # Required:
  #  - Output CSV report file option
  #
  # Next
  #  - Documentation
  #    + Main PDF/HTML, brief tutorial example
  #
  # To think about:
  #  - options for different RE strategies, e.g. no fill-in of sticky ends etc.
  #  - split input files to parallelise the initial stages of the process


if __name__ == '__main__':
  main_dir = os.path.dirname(os.path.dirname(__file__))
  sys.path.insert(0, main_dir)
  sys.path.insert(0, os.path.join(main_dir, 'nuc_tools')) 
  main()
