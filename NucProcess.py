import sys, os, re, shutil, gzip
import numpy as np

from collections import defaultdict
from subprocess import Popen, PIPE, call
from NucSvg import SvgDocument
from NucContactMap import nuc_contact_map


PROG_NAME = 'nuc_process'
VERSION = '1.0.0'
DESCRIPTION = 'Chromatin contact paired-read Hi-C processing module for Nuc3D and NucTools'
RE_CONF_FILE = 'enzymes.conf'
RE_SITES = {'MboI':'^GATC_', 'AluI':'AG^CT', 'BglII':'A^GATC_T',}
QUAL_SCHEMES = ['integer', 'phred33', 'phred64', 'solexa']
FASTQ_READ_CHUNK = 1048576
MIN_READ_LEN = 20
SCORE_TAG = re.compile(r'\sAS:i:(\S+)')
SCORE_TAG_SEARCH = SCORE_TAG.search
FILENAME_SPLIT_PATT = re.compile('[_\.]')
NCC_FORMAT = '%s %d %d %d %d %s %s %d %d %d %d %s %d %d %d\n'
  
  
if os.path.exists(RE_CONF_FILE):
  for line in open(RE_CONF_FILE):
    name, site = line.split()
    RE_SITES[name] = site


def open_file_r(file_path):

  if file_path.endswith('.gz'):
    file_obj = gzip.open(file_path, 'rt')
  else:
    file_obj = open(file_path, 'rU')

  return file_obj


def compress_file(file_path):

  in_file_obj = open(file_path)
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
   head_byte_1 = 0
   head_byte_2 = 0   
   
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
  

def remove_promiscuous(ncc_file, keep_files=True, zip_files=False, resolve_limit=1e3, close_cis=1e4):

  # resolve_limit : Allow two suitably close promiscous ends if long range cis or trans
 
  clean_ncc_file = tag_file_name(ncc_file, 'clean')
  
  if keep_files:
    promiscous_ncc_file = tag_file_name(ncc_file, 'promisc')  
    promisc_file_obj = open(promiscous_ncc_file, 'w')
    write_promisc = promisc_file_obj.write
  
  
  in_file_obj = open(ncc_file)
  clean_file_obj = open(clean_ncc_file, 'w')
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
  for re_start, re_ends in frag_counts.iteritems():
    if len(re_ends) > 1:
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
        
        else: # Other ends too far apart
          remove.add(re_start)
       
      else:
        remove.add(re_start)    
  
  n_trans = 0
  n_cis_near = 0
  n_cis_far = 0
      
  in_file_obj.seek(0)
  for line in in_file_obj:
    line_data = line.split()
    chr_a = line_data[0]
    chr_b = line_data[6]
    f_start_a = int(line_data[1])
    f_start_b = int(line_data[7])
    strand_a = line_data[5]
    strand_b = line_data[11]

    if (chr_a, f_start_a, strand_a) in remove:
      n_promiscuous += 1
      
      if keep_files:
        write_promisc(line)
      
    elif (chr_b, f_start_b, strand_b) in remove:
      n_promiscuous += 1
      
      if keep_files:
        write_promisc(line)
      
    else:
      n_clean += 1
      
      if chr_a != chr_b:
        n_trans += 1
      
      else:
        if strand_a == '+':
          p1 = int(line_data[2]) # Ligation junction is ahead at 3' end of RE frag
        else:
          p1 = f_start_a
          
        if strand_b == '+':
          p2 = int(line_data[8])
        else:
          p2 = f_start_b
     
        if abs(p1-p2) < 100000:
          n_cis_near += 1
 
        else:
          n_cis_far += 1
      
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
  
  promiscuity_stats = [('input_pairs', n),
                       ('promiscuous',(n_promiscuous, n)),
                       ('resolved',(n_resolved, n)),
                       ('clean',(n_clean, n))]

  final_stats = [('total_contacts', n),
                 ('cis_near',(n_cis_near, n)),
                 ('cis_far',(n_cis_far, n)),
                 ('trans',(n_trans, n))]


  return clean_ncc_file, promiscuity_stats, final_stats
  
 
def remove_redundancy(ncc_file, min_repeats=2, keep_files=True, zip_files=False, use_re_fragments=True):
  
  """
  A:B B:A redundancey taken care of at this stage because the NCC file pairs are internally sorted
  
  Option to remove redundancy at he RE fragment level (or otherwise at the read level) 
  
  Choose the longest read (not most common?) as the representitive for a merge
  """
   
  # First make temporary sorted file

  sort_file_name = tag_file_name(ncc_file, 'sort')
  cmd_args = ['sort', ncc_file]
  call(cmd_args, shell=False, stderr=None, stdin=None, stdout=open(sort_file_name, 'w'))
        
  sort_file_obj = open(sort_file_name, 'r')
  out_file_name = tag_file_name(ncc_file, 'multi_read')
  out_file_obj = open(out_file_name, 'w')
  
  if keep_files:
    uniq_file_name = tag_file_name(ncc_file, 'unique_read')
    uniq_file_obj = open(uniq_file_name, 'w')
  
  n_unique = 0
  n_redundant = 0
  n_pairs = 0
  mean_redundancy = 0.0
  
  # Compare using strand information given we want to know which fragment END is used
  if use_re_fragments:
    ncc_idx = (0,1,5,6,7,11) # Compare on RE fragment
  else:
    ncc_idx = (0,3,5,6,9,11) # Compare on read starts - does not consider ends as these can vary in read replicates
  
  line_prev = sort_file_obj.readline()
  line_data = line_prev.split()
  
  n = 1
  prev = [line_data[i] for i in ncc_idx]
  len_prev = abs(int(line_data[4]) - int(line_data[3])) + abs(int(line_data[10]) - int(line_data[9]))
    
  # Remove redundancy
  for line in sort_file_obj:
    n_pairs += 1
    line_data = line.split()
    curr = [line_data[i] for i in ncc_idx]
    len_curr = abs(int(line_data[4]) - int(line_data[3])) + abs(int(line_data[10]) - int(line_data[9]))
    
    if curr == prev:
      n += 1
      if len_curr > len_prev: # This repeat is better
        len_prev = len_curr
        line_prev = line
    
      continue
    
    else: # Not a repeat, write previous
      
      if n < min_repeats:
        if keep_files:
          uniq_file_obj.write(line_prev)
        n_unique += 1
      else:
        out_file_obj.write(line_prev)
        n_redundant += 1
      
      mean_redundancy += n
      
      line_prev = line
      len_prev = len_curr
      prev = curr
      n = 1
  
  # Write remaining    
  if n < min_repeats:
    if keep_files:
      uniq_file_obj.write(line_prev)
  else:
    out_file_obj.write(line_prev)
 
  out_file_obj.close()

  if keep_files:
    uniq_file_obj.close()
    
    if zip_files:
      compress_file(ncc_file)
  
  else:
    os.unlink(ncc_file)
  
  os.unlink(sort_file_name) # Remove temp file
  
  n = n_unique + n_redundant
  mean_redundancy /= float(n)
  
  stats = [('input_pairs', n_pairs),
           ('unique', (n_unique, n)),
           ('redundant', (n_redundant, n)),
           ('mean_redundancy', mean_redundancy)]
  
  return out_file_name, stats
  
  
def filter_pairs(pair_ncc_file, re1_file, re2_file=None, sizes=(100,2000), too_close_limit=1000, 
                 keep_files=True, zip_files=False, min_mappability=1, re2_tolerance=1000):
  
  re1_end_dict, re_end_dict1, chromo_name_dict1 = read_re_frag_file(re1_file)
  
  if re2_file:
    re2_frag_dict, re2_end_dict, chromo_name_dict = read_re_frag_file(re2_file)
  
  out_file_objs = {}
  out_file_names = {}
  counts = {}
  size_counts = defaultdict(int)
  min_size, max_size = sorted(sizes)
  frag_sizes = []
  
  for tag in ('accepted', 'too_small', 'too_big', 'circular_re1', 'internal_re1',
              'internal_re2', 'no_end_re2', 'overhang_re1', 'too_close',
              'adjacent_re1', 'low_mappability', 'unknown_contig'):
    counts[tag] = 0
    
    if keep_files or (tag == 'accepted'):
      file_name = tag_file_name(pair_ncc_file, 'filter_' + tag, '.ncc')
      file_obj = open(file_name, 'w')
      out_file_objs[tag] = file_obj
      out_file_names[tag] = file_name

  counts['near_cis_pairs'] = 0
  counts['far_cis_pairs'] = 0
  counts['trans_pairs'] = 0
  counts['input_pairs'] = 0
  
  def count_write(tag, line):
    counts[tag] += 1
    
    if keep_files or (tag == 'accepted'):
      out_file_objs[tag].write(line) 
  
  in_file_obj = open_file_r(pair_ncc_file)

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
    
    # Convert contig code names to proper chromosome names where possible
    chr_a = chromo_name_dict1.get(chr_a, chr_a)
    chr_b = chromo_name_dict1.get(chr_b, chr_b)
    
    if chr_a not in re_end_dict1:
      count_write('unknown_contig', line)
      continue  
 
    if chr_b not in re_end_dict1:
      count_write('unknown_contig', line)
      continue
    
    # Find which RE fragments the read positions wre within
    # index of frag with pos immediately less than end
    re1_a_idx = np.searchsorted(re_end_dict1[chr_a], [start_a])[0]
    re1_b_idx = np.searchsorted(re_end_dict1[chr_b], [start_b])[0]     
    
    re1_a_start, re1_a_end, mappability_a = re1_end_dict[chr_a][re1_a_idx] 
    re1_b_start, re1_b_end, mappability_b = re1_end_dict[chr_b][re1_b_idx]
    
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
    line = NCC_FORMAT % (chr_a, start_a, end_a, re1_a_start, re1_a_end, strand_a,
                         chr_b, start_b, end_b, re1_b_start, re1_b_end, strand_b,
                         int(ambig_group), int(pair_id), int(swap_pair))
 
    if (mappability_a < min_mappability) or (mappability_b < min_mappability):
      count_write('low_mappability', line)
      continue

    
    if re2_file:
      # Check Read is at the end of different RE2 fragments
      re2_a_idx = np.searchsorted(re2_end_dict[chr_a], [start_a])[0]
      re2_b_idx = np.searchsorted(re2_end_dict[chr_b], [start_b])[0]     
 
      # Internal RE2 not necessarily a problem if reads are in different RE1 fragements, e.g. potential ligation junction within
      if re2_a_idx == re2_b_idx and chr_a == chr_b and re1_a_idx == re1_b_idx:
        count_write('internal_re2', line)
        continue
 
      re2_a_start, re2_a_end, mappability_re2_a = re2_frag_dict[chr_a][re2_a_idx]
      re2_b_start, re2_b_end, mappability_re2_b = re2_frag_dict[chr_b][re2_b_idx]
       
      if (mappability_re2_a < min_mappability) or (mappability_re2_b < min_mappability):
        count_write('low_mappability', line)
        continue

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
        continue

      if abs(delta_re2_b) > re2_tolerance:
        count_write('no_end_re2', line)
        continue
   
    if chr_a == chr_b:
      
      if re1_a_idx == re1_b_idx: # Same Re1 fragment
        if start_a < start_b:
          if pos_strand_b and not pos_strand_a: # Reads go outwards in same frag: <-A B->
            count_write('circular_re1', line)
            continue
          
        else:
          if pos_strand_a and not pos_strand_b: # Reads go outwards in same frag: <-B A->
            count_write('circular_re1', line)
            continue      
    
        if (p1 <  re1_a_start < p2) or (p1 <  re1_a_end < p2): # Read A overshoots Re1 fragment end
          count_write('overhang_re1', line)

        elif (p3 < re1_b_start < p4) or (p3 < re1_b_end < p4): # Read B overshoots Re1 fragment end
          count_write('overhang_re1', line)
        
        else:
          count_write('internal_re1', line)

      elif abs(re1_a_idx-re1_b_idx) < 2: 
        count_write('adjacent_re1', line) # Mostly re-ligation
 
      else: # Different re1 fragment
        
        if pos_strand_a != pos_strand_b:
          
          if pos_strand_a and (p1 < p3): # Sequencing toward each other
            delta = p4 - p1 # separation of pair
            
            if delta < too_close_limit: 
              count_write('too_close', line)  # Pair is possible even without ligation 
              continue      
            
          elif pos_strand_b and (p3 < p1): # Sequencing toward each other
            delta = p2 - p3
          
            if delta < too_close_limit:
              count_write('too_close', line)  # Pair is possible even without ligation
              continue      

        if size_ok:
          if min(abs(p1-p3), abs(p2-p4)) < 1e4:
            counts['near_cis_pairs'] += 1
          else:
            counts['far_cis_pairs'] += 1

          count_write('accepted', line)
          
        elif size_t < min_size:
          count_write('too_small', line)
          
        else:
          count_write('too_big', line)
      
    else: # Good stuff, trans contact, many errors impossible
      if size_ok:
        counts['trans_pairs'] += 1
        count_write('accepted', line)
        
      elif size_t < min_size:
        count_write('too_small', line)
        
      else:
        count_write('too_big', line)
      
  for tag in out_file_objs:
    out_file_objs[tag].close()
    
  n = counts['input_pairs']
  sort_counts = []
  for key in ('input_pairs', 'accepted', 'near_cis_pairs', 'far_cis_pairs', 'trans_pairs',
              'internal_re1', 'adjacent_re1', 'circular_re1', 'overhang_re1',
              'too_close', 'too_small', 'too_big',
              'internal_re2', 'no_end_re2',
              'low_mappability', 'unknown_contig'):
   
    sort_counts.append((key, (counts[key], n))) # So percentages can be calculated in reports
  
  filter_file = out_file_names['accepted']  # Main output never compressed at this stage
  del out_file_names['accepted']
  
  if keep_files and zip_files:
    for tag in out_file_objs:
      compress_file(out_file_names[tag])    
          
  return filter_file, out_file_names, sort_counts, frag_sizes
    
    
def pair_sam_lines(line1, line2):
  """
  Convert from single to paired SAM format
  """
  
  read1 = line1.split('\t')
  read2 = line2.split('\t')
  
  #Relevant bitwise flags (flag in an 11-bit binary number)
  #1 The read is one of a pair
  #2 The alignment is one end of a proper paired-end alignment
  #4 The read has no reported alignments
  #8 The read is one of a pair and has no reported alignments
  #16 The alignment is to the reverse reference strand
  #32 The other mate in the paired-end alignment is aligned to the reverse reference strand
  #64 The read is the first (#1) mate in a pair
  #128 The read is the second (#2) mate in a pair

  #The reads were mapped as single-end data, so should expect flags of
  #0 (map to the '+' strand) or 16 (map to the '-' strand)
  #Output example: a paired-end read that aligns to the reverse strand
  #and is the first mate in the pair will have flag 83 (= 64 + 16 + 2 + 1)
  
  bits1 = int(read1[1])
  bits2 = int(read2[1])
  
  # Ignore non-alignments (and reads rejected due to -m 1 in Bowtie)
  if ( bits1 & 0x4 ) or ( bits2 & 0x4 ):
    return ( 0, 0 )
      
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
  
  
def pair_mapped_seqs(sam_file1, sam_file2, file_root, ambig=True, unique_map=False):
  
  paired_ncc_file_name = tag_file_name(file_root, 'pair', '.ncc')
  ambig_ncc_file_name = tag_file_name(file_root, 'pair_ambig', '.ncc')
  
  ncc_file_obj = open(paired_ncc_file_name, 'w')
  write_pair = ncc_file_obj.write

  ambig_ncc_file = open(ambig_ncc_file_name, 'w')
  write_ambig = ambig_ncc_file.write
    
  sort_sam_file1 = sam_file1 #_sort_sam_file_ids(sam_file1) # Only needed if Bowtie2 does not use --reorder option
  sort_sam_file2 = sam_file2 #_sort_sam_file_ids(sam_file2)
    
  file_obj1 = open_file_r(sort_sam_file1)
  file_obj2 = open_file_r(sort_sam_file2)
  
  n_map1  = 0
  n_map2  = 0
  n_pairs = 0
  n_ambig = 0
  n_unambig = 0
  n_unmapped = 0
     
  readline1 = file_obj1.readline
  readline2 = file_obj2.readline
          
  # Go through same files and pair based on matching id
  # Write out any to-many mapings to ambiguous
  
  def write_ncc(ncc_data_a, ncc_data_b, ambig_group, _id, ambig=False):
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
      
      
    if ambig:
      write_ambig(NCC_FORMAT % line_data)
    else:
      write_pair(NCC_FORMAT % line_data)
  
  id2 = 0
  line1 = readline1()
  line2 = readline2()  
  
  # Skip to end of headers
  
  while line1[0] == '@':
    line1 = readline1()

  while line2[0] == '@':
    line2 = readline2()
  
  # Process data lines
  
  id1 = line1[:10]
  id2 = line2[:10]
  
  while id1 and id2:
    _id = min(id1, id2)

    contact_a = []
    contact_b = []
    scores_a = []
    scores_b = []
    
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
    
    n_a = len(contact_a)
    n_b = len(contact_b)
    n = n_a * n_b
    n_pairs += 1
    
    if n < 1:
      n_unmapped += 1
    
    else:
      ambig = n > 1
       
      if n > 1:
        max_score_a = max(scores_a)
        max_score_b = max(scores_b)
        
        # Max score of zero means not less than perfect match
        if unique_map and max_score_a == 0 and max_score_b == 0 and scores_a.count(max_score_a) * scores_b.count(max_score_b) == 1:
          n_unambig += 1
          for ncc_a, score_a in contact_a:
            if score_a == max_score_a:
              for ncc_b, score_b in contact_b:
                if score_b == max_score_b:
                  write_ncc(ncc_a, ncc_b, n_pairs, int(_id), False)
                  break
              else:
                continue
              
              break  
              
        else:
          n_ambig += 1
          for ncc_a, score_a in contact_a:
            for ncc_b, score_b in contact_b:
              write_ncc(ncc_a, ncc_b, n_pairs, int(_id), True) # Write all pair combinations
            
      else:
        n_unambig += 1
        for ncc_a, score_a in contact_a:
          for ncc_b, score_b in contact_b:
            write_ncc(ncc_a, ncc_b, n_pairs, int(_id), False)
 
  #os.unlink(sort_sam_file1)
  #os.unlink(sort_sam_file2)
    
  stats = [('end_1_aligned', n_map1), 
           ('end_2_aligned', n_map2),
           ('unmapped_end', (n_unmapped, n_pairs)),
           ('unique',  (n_unambig, n_pairs)),
           ('ambiguous',    (n_ambig, n_pairs)),
           ('total_pairs',   n_pairs)]
           
  ncc_file_obj.close()
  
  return paired_ncc_file_name, ambig_ncc_file_name, stats


def map_reads(fastq_file, genome_index, align_exe, num_cpu, ambig=False, verbose=True):
  
  patt_1 = re.compile('(\d+) reads; of these:')
  patt_2 = re.compile('(\d+) \(.+\) aligned exactly 1 time')
  patt_3 = re.compile('(\d+) \(.+\) aligned 0 times')
  patt_4 = re.compile('(\d+) \(.+\) aligned >1 times')
  
  sam_file_path = tag_file_name(fastq_file, '', '.sam')
  
  
  cmd_args = [align_exe]
    
  cmd_args += [#'--very-sensitive',
               '-D', '20', '-R', '3', '-N', '0',  '-L', str(MIN_READ_LEN),  '-i', 'S,1,0.50',
               '-x', genome_index,
               '-k', '2',
               '--reorder',
               #'--un', fastq_file[:-6] + 'unaligned.fastq',  
               '-p', str(num_cpu),
               '-U', fastq_file,
               '-S', sam_file_path]
  
  n_reads = 0
  n_uniq  = 0
  n_unmap = 0
  n_ambig = 0
  
  proc = Popen(cmd_args, stdin=PIPE, stderr=PIPE)
  std_out, std_err = proc.communicate()
  
  if std_err:
    
    for line in std_err.split('\n'):
      
      #if verbose:
      print(line)
    
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
  
  return sam_file_path, stats
   

def clip_reads(fastq_file, file_root, junct_seq, replaced_seq, tag='clipped', min_len=MIN_READ_LEN):
  """
  Clips reads at ligation junctions, removes anbigous ends and discards very short reads
  """
  
  clipped_file = tag_file_name(file_root, tag, '.fastq')
  
  in_file_obj = open_file_r(fastq_file)
  n_rep = len(replaced_seq)
  n_junc = len(junct_seq)
  
  n_reads = 0
  n_clip = 0
  n_short = 0
  mean_len = 0
  
  out_file_obj = open(clipped_file, 'w', 2**16)
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
      n_clip += 1
      i = line2.index(junct_seq)
      line2 = line2[:i] + replaced_seq
      line4 = line4[:i+n_rep]

    while line2[-1] == 'N':
      line2 = line2[:-1]
      line4 = line4[:-1]
    
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
    mean_len /= n_reads

  stats = [('input_reads',n_reads),
           ('clipped',(n_clip, n_reads)),
           ('too_short',(n_short, n_reads)),
           ('mean_length',mean_len)]
  
  return clipped_file, stats
  

def get_chromo_re_fragments(sequence, re_site, offset, genome_index, align_exe, num_cpu, mappability_length=50):
  
  frag_data = []
  frag_data_append = frag_data.append
  re_site = re_site.replace('^', '')
  re_site = re_site.replace('_', '')
  
  if 'N' in re_site:
    re_site = re_site.replace('N', '[AGCT]') # Could allow more ambiguiuty codes, if feeling generous
    frags = re.split(re_site, sequence)
  else:
    frags = sequence.split(re_site)
  
  # Make FASTA file of fragments to assess mappability
  fasta_file_name = genome_index + '_temp_re_frags.fasta'
  fasta_file = open(fasta_file_name, 'w')
  fasta_write = fasta_file.write
  step = mappability_length/2
  site_len = len(re_site)
  offset_start = site_len - offset

  frag_start = offset_start
  for frag_seq in frags:
    n = len(frag_seq)
    
    start = frag_start - offset_start
    end   = frag_start + n + offset # Offset is cut site, which can be in subseq removed during split
    
    if not n:
      frag_start += site_len
      continue
     
    a = 0
    b = mappability_length
    
    while b < n:
      fasta_line = '>%s_%s\n%s\n' % (start, a, frag_seq[a:b])
      fasta_write(fasta_line)
      a += step
      b = a + mappability_length
    
    if b != n:
      a = n-mappability_length
      b = n
      fasta_line = '>%s_%s\n%s\n' % (start, a, frag_seq[a:b])
      fasta_write(fasta_line)
    
    frag_data_append([start, end-1]) # Mappability added later
    frag_start += n + site_len

  frag_data[-1][1] = len(sequence)-1 # Last position is seq end, not cut site
  
  fasta_file.close()
  
  # Calculate mappability
  
  frag_read_counts = defaultdict(int)
  sam_file_path = tag_file_name(fasta_file_name, '', '.sam')

  cmd_args = [align_exe, '--sensitive',
              '-x', genome_index,
              '-p', str(num_cpu),
              '-f', '-U', fasta_file_name,
              '-S', sam_file_path]  # Unclude unaligned fragments, naturally
  
  #print cmd_args
  
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
    start, indent = fasta_head.split('_')

    frag_read_counts[int(start)] += 1
    
  sam_file_obj.close()
  os.unlink(sam_file_path)
  os.unlink(fasta_file_name)
  
  # Add mappability
  for i, (start, end) in enumerate(frag_data):
    frag_data[i] = [start, end, frag_read_counts[start]]
  
  return frag_data
    

def get_re_seq(re_name):

  site = RE_SITES[re_name]
  return site.replace('^', '').replace('_','')
  

def check_re_frag_file(genome_index, re_name, genome_fastas, align_exe, num_cpu, verbose=False, remap=False):
  
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
  
  site = RE_SITES[re_name]
  seq = get_re_seq(re_name)
  pos = site.replace('_', '').index('^')
  
  out_file_obj = open(re_site_file, 'w')
  
  head = 'Chromo\tContig\tFrag_Start\tFrag_End\tMappability\n'
  msg = 'Calculating %s fragment locations and mappability for chromosome %s contig %s'
  out_file_obj.write(head)
  
  for fasta_file in genome_fastas:
    file_obj = open(fasta_file)
    
    #chromo = None
    contig = None
    chromo = os.path.splitext(fasta_file)[0].split('_')[-1]
    seq = ''
    
    for fasta_line in file_obj:
      if fasta_line[0] == '>':
        if chromo and seq: # add previous
          if verbose:
            info(msg % (re_name, chromo, contig))
        
          for start, end, mappability in get_chromo_re_fragments(seq, site, pos, genome_index, align_exe, num_cpu):
            frag_line = '%s\t%s\t%d\t%d\t%d\n' % (chromo, contig, start, end, mappability)
            out_file_obj.write(frag_line)
            
        contig = fasta_line[1:].split()[0]
        seq = ''
      
      else:
        seq += fasta_line[:-1]
    
    file_obj.close()
    
    if chromo and seq:
      if verbose:
        info(msg % (re_name, chromo, contig))
        
      for start, end, mappability in get_chromo_re_fragments(seq, site, pos, genome_index, align_exe, num_cpu):
        frag_line = '%s\t%s\t%d\t%d\t%d\n' % (chromo, contig, start, end, mappability)
        out_file_obj.write(frag_line)
 
       
  out_file_obj.close()
  
  return re_site_file
  

def read_re_frag_file(file_path):

  file_obj = open_file_r(file_path)
  file_obj.readline()
  
  data = [line.split() for line in file_obj.readlines()]
  
  frag_dict = defaultdict(list)
  end_dict = defaultdict(list)
  chromo_contigs = defaultdict(set)
  
  for chromo, contig, start, end, mappability in data:
    start = int(start)
    end = int(end)
    mappability = int(mappability)
    frag_dict[contig].append( (start, end, mappability) )  
    end_dict[contig].append(end)
    chromo_contigs[chromo].add(contig)
  
  chromo_name_dict = {} # Safely mappable contig names
  for chromo in chromo_contigs:
    contigs = chromo_contigs[chromo]
    
    if len(contigs) == 1:
      contig = contigs.pop()
      chromo_name_dict[contig] = chromo
      chromo_name_dict[chromo] = chromo
      
      frag_dict[chromo] = frag_dict[contig]
      end_dict[chromo] = end_dict[contig] = np.array(end_dict[contig])
    
    else:
      for contig in contigs:
        end_dict[contig] = np.array(end_dict[contig])
      

  return dict(frag_dict), dict(end_dict), chromo_name_dict
  

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

  fasta_file_str= ','.join(fasta_files)

  cmd_args = [indexer_exe, '-f']
  
  if quiet:
    cmd_args.append('-q')
  
  if pack:
    cmd_args.append('-p')

  cmd_args += ['-t', str(table_size), fasta_file_str, base_name]

  #print(' '.join(cmd_args))
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


def is_genome_indexed(file_path, file_ext='.1.bt2'):
 
  return os.path.exists(file_path + file_ext)
  

def check_index_file(file_path, file_exts=('.1.bt2', '.2.bt2', '.3.bt2', '.4.bt2', '.rev.1.bt2', '.rev.2.bt2')):
  
  for file_ext in file_exts:
    is_ok, msg = check_regular_file(file_path + file_ext)
    
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
    msg = 'File "%s" is of zero size '% file_path
    return False, msg
    
  if not os.access(file_path, os.R_OK):
    msg = 'File "%s" is not readable' % file_path
    return False, msg
  
  return True, msg
  
  
def info(msg, prefix='INFO'):

  print('%8s : %s' % (prefix, msg))


def warn(msg, prefix='WARNING'):

  print('%8s : %s' % (prefix, msg))


def fatal(msg, prefix='%s FAILURE' % PROG_NAME):

  print('%8s : %s' % (prefix, msg))
  print('Use %s with -h to display command line options' % PROG_NAME)
  sys.exit(0)


def testImports():
  
  critical = False
  
  try:
    import samtools
    
  except ImportError:
    critical = True
    warn('pysam module not available')

  if critical:
    sys.exit(0)


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
    scheme = 'integer-quals'
  
  elif (max_qual < 75) and (min_qual < 59): # Sanger, Illumina 1.8+
    scheme = 'phred33-quals'
  
  elif (max_qual > 74): 
    if min_qual < 64:
      scheme = 'solexa-quals'
    else:
      scheme = 'phred64-quals'
  
  else:
    warn('FASTQ quality scheme could not be determined. Assuming Phred+33 (Illumina 1.8+)')
    scheme = 'phred33-quals'

  return scheme

def write_report(report_file, frag_sizes, general_stats, clip_stats1, clip_stats2, sam_stats1, sam_stats2,
                 pair_stats, filter_stats, redundancy_stats, promiscuity_stats, final_stats):
  
  def format_list(d):
    
    l = []
    for k, v in d:
      c1 = k.replace('_', ' ')
      c1 = c1[0].upper() + c1[1:]
           
      if isinstance(v, (tuple, list)):
        v, n = v
        percent = 100.0 * v/float(n)
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
    
    return sizes, names_2
      
      
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
  
  y += head_height
  svg_doc.text('Pairing reads', (x,y), head_size, bold=True, color=head_color)

  y += head_pad
  data = format_list(pair_stats)
  tw, th = svg_doc.table(x, y, table_width/2, data, False, text_anchors, col_formats=None, size=row_size, main_color=main_color)
  
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
           'internal_re2','no_end_re2', 'low_mappability', 'unknown_contig']
  colors = ['#80C0FF','#FFFF00','#FF0000','#FFC000','#FF0080',
            '#00FFFF','#00C0C0','#008080',
            '#40FF40','#008000','#FF00FF','#404040']
  vals, names = pie_values(filter_stats, names)
  svg_doc.pie_chart(x1, y, chart_height, vals, names, colors, line_color='#808080', small_val=0.01)
  
  hist, edges = np.histogram(frag_sizes, bins=50, range=(0, 1100))
  data_set = [(int(edges[i]), int(val)) for i, val in enumerate(hist)]
  svg_doc.graph(x1, y+chart_height, table_width/2, th-chart_height-40, [data_set], 'Size', 'Count',
                names=None, colors=None,  graph_type='line',
                symbols=None, line_widths=None, symbol_sizes=None,
                legend=False, title=None, plot_offset=(80, 20))
  y += th
  
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
  names = ['trans','cis_far','cis_near']
  colors = ['#80C0FF','#FFFF00','#FF0000']
  vals, names = pie_values(final_stats, names)
  svg_doc.pie_chart(x1, y, chart_height, vals, names, colors, line_color='#808080')
  y += chart_height
  
  
  height = y + main_pad
  width = table_width + 3 * main_pad
  
  svg_doc.write_file(report_file, width, height)



def log_report(report_file, section, data, verbose=False):
  
  def format_val(val):
    if isinstance(val, int):
      v = '{:,}'.format(val)
    
    elif isinstance(val, float):
      v = '{:.3f}'.format(val)
   
    else:
      v = '{}'.format(val)
    
    return v
  
  lines = [section]

  for key, val in data:
    name = key.replace('_',' ')
    
    if isinstance(val, (tuple, list)):
      val, n = val
      percent = 100.0 * val/float(n)
      info_line = '  {} : {} ({:.2f}%)'.format(key, val, percent)
    
    else:
      val = format_val(val)
      info_line = '  {} : {}'.format(key, val)

    lines.append(info_line)

  if verbose:
    for line in lines:
      info(line)

def nuc_process(fastq_paths, genome_index, re1, re2=None, sizes=(300,800), min_rep=2, num_cpu=1,
                ambig=True, unique_map=False, out_file=None, ambig_file=None, report_file=None,
                align_exe=None, qual=None, g_fastas=None, is_pop_data=False, remap=False, reindex=False,
                keep_files=True, lig_junc=None, zip_files=True, sam_format=True, verbose=True):
  """
  Main function for command-line operation
  """
  
  # Files can be empty if just re-indexing etc...
  if not fastq_paths and not (reindex and g_fastas):
    fatal('No FASTQ files specified')
  
  # Ambiguous data is only for single-cell use
  if is_pop_data and ambig:
    fatal('ambiguous option is incompatible with population Hi-C data')
  
  # Check FASTQ files : two present, exist, right format
  
  if fastq_paths:
    if len(fastq_paths) == 1:
      fatal('Only one FASTQ file specified (exactly two required)')

    if len(fastq_paths) > 2:
      fatal('More than two FASTQ files specified (exactly two required)')

  for file_path in fastq_paths:
    is_ok, msg = check_regular_file(file_path)
    
    if not is_ok:
      fatal(msg)
    
    is_ok, msg = check_fastq_file(file_path)
    
    if not is_ok:
      fatal(msg)

  # Check genome index file, if not being rebuilt
  if is_genome_indexed(genome_index):
    is_ok, msg = check_index_file(genome_index)
    
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
   
  # Check restriction enzymes  
  if re1 not in RE_SITES:
    msg = 'Restriction enzyme "%s" not known. Available: %s.' % (re1, ', '.join(sorted(RE_SITES)))
    msg += ' %s can be edited to add further cut-site defintions.' % RE_CONF_FILE
    fatal(msg)
    
  is_ok, msg = check_re_site(RE_SITES[re1])
  if not is_ok:
    fatal(msg)
  
  re1Seq = get_re_seq(re1)
  
  if re2:
    re2Seq = get_re_seq(re2)
    
    if re2 not in RE_SITES:
      msg = 'Restriction enzyme "%s" not known. Available: %s.' % (re2, ', '.join(sorted(RE_SITES)))
      msg += ' Enzymes.conf can be edited to add further cut-site defintions.'
      fatal(msg)
  
    is_ok, msg = check_re_site(RE_SITES[re2])
    if not is_ok:
      fatal(msg)
      
  else:
    re2Seq = None  
    
  # Check FASTQ quality scheme
  if qual:
    if qual not in QUAL_SCHEMES:
      msg = 'FASTQ quality scheme "%s" not known. Available: %s.' % (qual, ', '.join(sorted(QUAL_SCHEMES)))
      msg += ' Scheme will be deduced automatically if not specified.'
      fatal(msg)
    
    qual = qual + '_quals'
  
  elif fastq_paths:
    qual = get_fastq_qual_scheme(fastq_paths[0])
    
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
    file_roots = []
    for fastq_path in fastq_paths:
      if fastq_path.lower().endswith('.gz'):
        fastq_path = fastq_path[:-3]
    
      dir_name, file_name = os.path.split(fastq_path)
      file_roots.append(os.path.splitext(file_name)[0])
 
    file_root = '_'.join(file_roots)
    file_root = os.path.join(dir_name, file_root)
  
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
    contact_map_file = os.path.splitext(report_file)[0] + '_contact_map.svg'
  else:
    report_file = file_root + '_report.svg'
    contact_map_file =  file_root + '_contact_map.svg'

  # Check for no upper fragment size limit
  if len(sizes) == 1:
    size_str = '%d - unlimited' % tuple(sizes)
    sizes = (sizes[0], int(1e16))
  else:
    size_str = '%d - %d' % tuple(sizes)
  
  general_stats = [('Input FASTQ file 1',fastq_paths[0]),
                   ('Input FASTQ file 2',fastq_paths[1]),
                   ('Output file (main)',os.path.abspath(out_file)),
                   ('Output file (ambiguous)',os.path.abspath(ambig_file) if ambig_file else 'None'),
                   ('Report file',os.path.abspath(report_file)),
                   ('Intermediate file directory',os.path.abspath(intermed_dir)),
                   ('Aligner executable',align_exe),
                   ('Genome index',genome_index),
                   ('FASTQ quality scheme',qual),
                   ('RE1 site',re1),
                   ('RE2 site',(re2 or 'None')),
                   ('Ligation junction',lig_junc),
                   ('Molecule size range',size_str),
                   ('Min sequencing repeats',min_rep),
                   ('Parallel CPU cores',num_cpu),
                   ('Keep intermediate files?', 'Yes' if keep_files else 'No'),
                   ('Input is single-cell Hi-C?', 'No' if is_pop_data else 'Yes'),
                   ('SAM output?', 'Yes' if sam_format else 'No'),
                   ('GZIP output?', 'Yes' if zip_files else 'No'),
                  ]

  log_report(report_file, 'General', general_stats, verbose)
  
  # Check if genome needs indexing  
  if g_fastas and (reindex or not is_genome_indexed(genome_index)):
    indexer_exe = os.path.join(os.path.dirname(align_exe), 'bowtie2-build')
    
    warn('Indexing genome, this may take some time...')
    output_dir, base_name = os.path.split(genome_index)
    index_genome(base_name, g_fastas, output_dir, indexer_exe) # Latest version of bowtie2 can do parallel index builds (--threads)
  
  
  # Create RE fragments file if not present
  re1_file = check_re_frag_file(genome_index, re1, g_fastas, align_exe, num_cpu, verbose=verbose, remap=remap)
  
  if re2:
    re2_file = check_re_frag_file(genome_index, re2, g_fastas, align_exe, num_cpu, verbose=verbose, remap=remap)
  else:
    re2_file = None

  # Clip read seqs at any sequenced ligation junctions
  if verbose:
    info('Clipping FASTQ reads...')

  clipped_file1, clip_stats1 = clip_reads(fastq_paths[0], intermed_file_root, lig_junc, replaced_seq=re1Seq, tag='reads1_clipped')
  log_report(report_file, 'Clip Reads 1', clip_stats1, verbose)
  
  clipped_file2, clip_stats2 = clip_reads(fastq_paths[1], intermed_file_root, lig_junc, replaced_seq=re1Seq, tag='reads2_clipped')
  log_report(report_file, 'Clip Reads 2', clip_stats2, verbose)
 
  # Do the main genome mapping
  if verbose:
    info('Mapping FASTQ reads...')
 
  sam_file1, sam_stats1 = map_reads(clipped_file1, genome_index, align_exe, num_cpu, ambig, verbose=verbose)
  log_report(report_file, 'Map Reads 1', sam_stats1, verbose)
 
  sam_file2, sam_stats2 = map_reads(clipped_file2, genome_index, align_exe, num_cpu, ambig, verbose=verbose)
  log_report(report_file, 'Map Reads 2', sam_stats2, verbose)
   
  if verbose:
    info('Pairing FASTQ reads...')
    
  paired_ncc_file, ambig_paired_ncc_file, pair_stats = pair_mapped_seqs(sam_file1, sam_file2, intermed_file_root, ambig, unique_map)
  log_report(report_file, 'Read pairing', pair_stats, verbose)
  
  # Write SAM, if requested
  if sam_format:
    write_sam_file(paired_ncc_file, sam_file1, sam_file2)
    
    if ambig:
      write_sam_file(ambig_paired_ncc_file, sam_file1, sam_file2)
  
  # Filtering
  if verbose:
    info('Filtering mapped sequences...')
  
  filter_output = filter_pairs(paired_ncc_file, re1_file, re2_file, sizes, keep_files, zip_files)
  filter_ncc_file, fail_file_names, filter_stats, frag_sizes = filter_output
  
  log_report(report_file, 'Fragment filtering', filter_stats, verbose)
  
  # Write SAM, if requested
  if sam_format:
    write_sam_file(filter_ncc_file, sam_file1, sam_file2)
    
    for key in fail_file_names:
      write_sam_file(fail_file_names[key], sam_file1, sam_file2)
      
  if ambig:
    if verbose:
      info('Filtering ambiguously mapped sequences...')
  
    filter_output = filter_pairs(ambig_paired_ncc_file, re1_file, re2_file, sizes, keep_files, zip_files)
    ambig_filter_ncc_file, ambig_fail_file_names, ambig_filter_stats, ambig_frag_sizes = filter_output
    
    log_report(report_file, 'Fragment filtering (ambiguous)', ambig_filter_stats, verbose)
 
    # Write SAM, if requested
    if sam_format:
      write_sam_file(ambig_filter_ncc_file, sam_file1, sam_file2)
 
      for key in fail_file_names:
        write_sam_file(ambig_fail_file_names[key], sam_file1, sam_file2)
      
  # Merge duplicates
  if verbose:
    info('Removing duplicate contacts...')
  
  nr_ncc_file, redundancy_stats = remove_redundancy(filter_ncc_file, min_rep, keep_files, zip_files)
  log_report(report_file, 'Duplicate removal', redundancy_stats, verbose)
  
  if sam_format:
    write_sam_file(nr_ncc_file, sam_file1, sam_file2)

  if ambig:
    if verbose:
      info('Filtering duplicate ambiguous contacts...')
      
    ambig_nr_ncc_file, ar_stats = remove_redundancy(ambig_filter_ncc_file, min_rep, keep_files, zip_files)
    log_report(report_file, 'Duplicate removal (ambiguous)', ar_stats, verbose)
 
    if sam_format:
      write_sam_file(ambig_nr_ncc_file, sam_file1, sam_file2)
    
  # Remove promiscuous ends
  
  if is_pop_data:
    promiscuity_stats = {}
    shutil.copyfile(nr_ncc_file, out_file)
    
    if sam_format:
      write_sam_file(nr_ncc_file, sam_file1, sam_file2)
    
    if keep_files:
      if zip_files:
        nr_ncc_file = compress_file(nr_ncc_file)
      
    else:
      os.unlink(nr_ncc_file)
 
  else:
    if verbose:
      info('Removing pairs with promiscuous ends...')
 
    clean_ncc_file, promiscuity_stats, final_stats = remove_promiscuous(nr_ncc_file, keep_files, zip_files)
    log_report(report_file, 'Promiscuous end removal', promiscuity_stats, verbose)

    if sam_format:
      write_sam_file(clean_ncc_file, sam_file1, sam_file2)
 
    shutil.copyfile(clean_ncc_file, out_file)
 
    if ambig:
      if verbose:
        info('Removing ambiguous pairs with promiscuous ends...')
 
      ambig_clean_ncc_file, ap_stats, afs = remove_promiscuous(ambig_nr_ncc_file, keep_files, zip_files)
      log_report(report_file, 'Promiscuous end removal (ambiguous)', ap_stats, verbose)

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
  
  write_report(report_file, frag_sizes, general_stats, clip_stats1, clip_stats2, sam_stats1, sam_stats2,
               pair_stats, filter_stats, redundancy_stats, promiscuity_stats, final_stats)

  nuc_contact_map(out_file, contact_map_file)

  if verbose:
    info('Nuc Process all done.')

  return out_file
  


if __name__ == '__main__':
 
  from argparse import ArgumentParser
  
  epilog = 'Note %s can be edited to add further restriction enzyme cut-site defintions. ' % RE_CONF_FILE
  epilog += 'For further help email tjs23@cam.ac.uk or wb104@cam.ac.uk'
  
  res = sorted(RE_SITES)
  avail_re = 'Available: ' + ', '.join(res)
  avail_quals = 'Available: ' + ', '.join(QUAL_SCHEMES)
  
  arg_parse = ArgumentParser(prog=PROG_NAME, description=DESCRIPTION,
                            epilog=epilog, prefix_chars='-', add_help=True)

  arg_parse.add_argument('-i', nargs=2, metavar='FASTQ_FILE',
                         help='Input paired-read FASTQ files to process (accepts wildcards that match two files)')

  arg_parse.add_argument('-g',  metavar='GENOME_FILE',
                         help='Genome index file to map sequence reads to. A new index will be created with this name if the index is missing and genome FASTA files are specified')

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
                         help='Optional output name for NCC format chromosome contact file')

  arg_parse.add_argument('-oa', metavar='NCC_FILE',
                         help='Optional output name for ambiguous contact NCC file')

  arg_parse.add_argument('-or', metavar='REPORT_FILE',
                         help='Optional output name for SVG format report file')

  arg_parse.add_argument('-b', metavar='EXE_FILE',
                         help='Path to bowtie2 (read aligner) executable (will be searched for if not specified)')

  arg_parse.add_argument('-q', metavar='SCHEME',
                         help='Use a specific FASTQ quality scheme (normally not set and deduced automatically). ' + avail_quals)

  arg_parse.add_argument('-m', default=False, action='store_true',
                         help='Force a re-mapping of genome restriction enzyme sites (otherwise cached values will be used if present)')

  arg_parse.add_argument('-p', default=False, action='store_true',
                         help='The input data is population Hi-C; single-cell processing steps are avoided')

  arg_parse.add_argument('-x', '--reindex', default=False, action='store_true', dest='x',
                         help='Force a re-indexing of the genome (given appropriate FASTA files)')

  arg_parse.add_argument('-f', nargs='+', metavar='FASTA_FILES',
                         help='Specify genome FASTA files for index building (accepts wildcards)')

  arg_parse.add_argument('-a', default=False, action='store_true',
                         help='Whether to report ambiguously mapped contacts')
                        
  arg_parse.add_argument('-k', default=False, action='store_true',
                         help='Keep any intermediate files (e.g. clipped FASTQ etc). Note initial, primary SAM files are always kept.')

  arg_parse.add_argument('-sam', default=False, action='store_true',
                         help='Write paired contacts files to SAM format')

  arg_parse.add_argument('-l', metavar='SEQUENCE',
                         help='Seek a specific ligation junction sequence (otherwise this is guessed from the primary restriction enzyme)')

  arg_parse.add_argument('-z', default=False, action='store_true',
                         help='GZIP compress any output FASTQ files')
                        
  arg_parse.add_argument('-v', '--verbose', action='store_true', dest='v',
                         help='Display verbose messages to report progress')    
 
  arg_parse.add_argument('-u', default=False, action='store_true',
                         help='Whether to only accept uniquely mapping genome positions and not attempt to resolve certain classes of ambiguous mapping.')
  
  args = vars(arg_parse.parse_args())
  
  fastq_paths  = args['i']
  genome_index = args['g']
  re1          = args['re1']
  re2          = args['re2']
  sizes        = args['s']
  num_cpu      = args['n']
  min_rep      = args['r']
  ambig        = args['a']
  out_file     = args['o']
  ambig_file   = args['oa']
  report_file  = args['or']
  align_exe    = args['b']
  qual         = args['q']
  g_fastas     = args['f']
  is_pop_data  = args['p']
  remap        = args['m']
  reindex      = args['x']
  keep_files   = args['k']
  lig_junc     = args['l']
  zip_files    = args['z']
  unique_map   = args['u']
  sam_format   = args['sam']
  verbose      = args['v']
  
  if sizes:
    sizes = sorted([int(x) for x in re.split('\D+', sizes)])
  
  nuc_process(fastq_paths, genome_index, re1, re2, sizes, min_rep, num_cpu,
              ambig, unique_map, out_file, ambig_file, report_file, align_exe,
              qual, g_fastas, is_pop_data, remap, reindex, keep_files,
              lig_junc, zip_files, sam_format, verbose)
    
  # Required:
  #  - hist on-the-fly clipped read and frag size length distributions
  #  - Output CSV report file option
  #
  # Next
  #  - Documentation
  #    + Main PDF/HTML
  #    + Tutorial
  #    + Example data
  #
  # To think about:
  #  - could clip sequences on a quality score threshold
  #  - could add an option to separate isolated contacts
  #  - barcode demultiplexing should be a related but separate program
  #  - options for different RE strategies, e.g. no fill-in of sticky ends etc.
  #  - population option not exclude unsupported/single read pairs,  low-mappabaility regions
  #  - normalise/correct population data?
