"""
---- COPYRIGHT ----------------------------------------------------------------

Copyright (C) 20016-2026
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
import os, sys, shutil
from math import log, exp, ceil
from time import time
from collections import defaultdict

import numpy as np
from scipy import signal

from hic_core.nuc_process import info, warn, fatal

PROG_NAME = 'single_cell_hic_disambiguate'
VERSION = '1.1.0'
DESCRIPTION = 'An extension tool for nuc_process that attempts to resolve ambiguous (multi-mapping) contacts in single-cell Hi-C data. Operates on NCC format files.'

MIN_SCORE_SEQ_SEP = 100
ISOLATION_THRESHOLD = 0.01

DEFAULT_SEP_THRESHOLD  = 4e7
DEFAULT_SCORE_THRESHOLD = 10.0
DEFAULT_MIN_CONTACTS = int(1e5)
DEFAULT_OUT_SUFFIX = '_disambig'
DEFAULT_SMALL_CHR_SIZE_LIMIT = 1.1e8 # Larger than this can have isolated unambiguous (non-homologue) trans contacts removed


def _get_network_score(sep_threshold, chr_a, chr_b, pos_a, pos_b, bin_a, bin_b,
                       unambig_contacts, primary_limit=8.0, min_sep=MIN_SCORE_SEQ_SEP):
  """
  primary_limit is an early-stopping heuristic to improve performance; should not be changed
  """
  
  sep_scale = 0.5 * sep_threshold
  plim = 10.0 ** primary_limit
  score_l = 1.0
  score_u = 1.0
  unambig_contacts_get = unambig_contacts.get
  
  a_bins = (bin_a-1, bin_a, bin_a+1)
  b_bins = (bin_b-1, bin_b, bin_b+1)
  
  for a in a_bins:
    rev = a < bin_a
    peri_a = a != bin_a
    
    for b in b_bins:
      key2 = (chr_a, chr_b, a, b)
      
      unambig = unambig_contacts_get(key2)
      
      if unambig:
        if rev:
          unambig = unambig[::-1]
        
        if bin_a == a and bin_b == b:
          for pos_1, pos_2 in unambig:
            delta_1 = pos_1-pos_a
            delta_2 = pos_2-pos_b
 
            if delta_1 and delta_2:
              if delta_1 < 0.0:
                delta_1 = max(-delta_1, min_sep)/sep_scale
                score_l += exp(-delta_1*delta_1)
              else:
                delta_1 = max(delta_1, min_sep)/sep_scale
                score_u += exp(-delta_1*delta_1)

              if delta_2 < 0.0:
                delta_2 = max(-delta_2, min_sep)/sep_scale
                score_l += exp(-delta_2*delta_2)
              else:
                delta_2 = max(delta_2, min_sep)/sep_scale
                score_u += exp(-delta_2*delta_2)
         
        else:
          for pos_1, pos_2 in unambig:
            delta_1 = pos_1-pos_a
            delta_2 = pos_2-pos_b
 
            if (abs(delta_1) < sep_threshold) and (abs(delta_2) < sep_threshold):
 
              if delta_1 < 0.0:
                delta_1 = max(-delta_1, min_sep)/sep_scale
                score_l += exp(-delta_1*delta_1)
              else:
                delta_1 = max(delta_1, min_sep)/sep_scale
                score_u += exp(-delta_1*delta_1)

              if delta_2 < 0.0:
                delta_2 = max(-delta_2, min_sep)/sep_scale
                score_l += exp(-delta_2*delta_2)
              else:
                delta_2 = max(delta_2, min_sep)/sep_scale
                score_u += exp(-delta_2*delta_2)
 
            elif peri_a and abs(delta_1) > sep_threshold:
              break
        
        if score_l * score_u > plim:
          return primary_limit

  return log(score_l * score_u, 10.0)


def _write_ambig_filtered_ncc(in_file_path, out_ncc_path, ag_data, resolved_ag=None, removed_ag=None):
  
  if not resolved_ag:
    resolved_ag = {}
    
  if not removed_ag:
    removed_ag = set()
  
  n = 0
  ambig_group = 0
  seen_res_group = set()
  
  with open(out_ncc_path, 'w') as out_file_obj, open(in_file_path) as in_file_obj:
    write = out_file_obj.write

    for i, line in enumerate(in_file_obj):
      row = line.split()
      chr_a = row[0]
      chr_b = row[6]
      ambig_code = row[12]

      if chr_a > chr_b:
        chr_a, chr_b = chr_b, chr_a

      if '.' in ambig_code: # Updated NCC format
        sz = int(float(ambig_code))
        if sz > 0:
          ambig_group += 1
          
      else:
        ambig_group = int(ambig_code)    
        sz = len(ag_data[ambig_group])
     
      if ambig_group in resolved_ag:    
        keep = 1 if i in resolved_ag[ambig_group] else 0 # Else inactivate
        
        if ambig_group in seen_res_group:
          row[12] = f'0.{keep}'
        else:
          seen_res_group.add(ambig_group)
          row[12] = f'{sz}.{keep}'

      elif ambig_group in removed_ag: # Inactivate
        if ambig_group in seen_res_group:
          row[12] = '0.0'
        else:
          seen_res_group.add(ambig_group)
          row[12] = f'{sz}.0'
           
      line = ' '.join(row)+'\n'
      write(line)    
      
      n += 1

  msg = " .. written {:,} of {:,} lines to {}".format(n, i+1, out_ncc_path)
  
  info(msg)
  
      
def _load_bin_sort_ncc(in_ncc_path, sep_threshold):
  
  msg = ' .. reading contacts from %s' % in_ncc_path
  info(msg)

  ag_data = defaultdict(list)
  nonpa_contacts = defaultdict(list)
  ag_positions = defaultdict(set)
  pos_ags = set()
  
  # Identify positional ambiguity
  
  with open(in_ncc_path) as in_file_obj:
    ambig_group = 0

    for line in in_file_obj:
      chr_a, start_a, end_a, f_start_a, f_end_a, strand_a, chr_b, start_b, end_b, \
        f_start_b, f_end_b, strand_b, ambig_code, pair_id, swap_pair = line.split()
        
      if '.' in ambig_code: # Updated NCC format
        if int(float(ambig_code)) > 0:
          ambig_group += 1 # Count even if inactive; keep consistent group numbering

      else:
        ambig_group = int(ambig_code)          

      if ambig_code.endswith('.0'): # Inactive
        continue
      
      if strand_a == '+':
        pos_a = int(f_end_a)
      else:
        pos_a = int(f_start_a)

      if strand_b == '+':
        pos_b = int(f_end_b)
      else:
        pos_b = int(f_start_b)
      
      ag_positions[ambig_group].add((chr_a, pos_a))
      ag_positions[ambig_group].add((chr_b, pos_b))
  
  for ag in ag_positions:
    chr_pos = ag_positions[ag]
    
    if len(chr_pos) > 4:
      pos_ags.add(ag)
      
    else:
      chr_roots = set([x[0].split('.')[0] for x in chr_pos])
      
      if len(chr_roots) > 2:
        pos_ags.add(ag)
  
  msg = ' .. found {:,} positional ambiguity groups from {:,} total'.format(len(pos_ags), ambig_group)
  info(msg)
    
  # Bin and group
  pa_contacts = defaultdict(list)
  i = -1 # File could be empty
  
  with open(in_ncc_path) as in_file_obj:
    ambig_group = 0

    for i, line in enumerate(in_file_obj):
      chr_a, start_a, end_a, f_start_a, f_end_a, strand_a, chr_b, start_b, end_b, \
        f_start_b, f_end_b, strand_b, ambig_code, pair_id, swap_pair = line.split()

      if '.' in ambig_code: # Updated NCC format
        if int(float(ambig_code)) > 0:
          ambig_group += 1
      else:
        ambig_group = int(ambig_code)          
      
      if ambig_code.endswith('.0'): # Inactive
        continue
      
      if strand_a == '+':
        pos_a = int(f_end_a)
      else:
        pos_a = int(f_start_a)

      if strand_b == '+':
        pos_b = int(f_end_b)
      else:
        pos_b = int(f_start_b)
 
      if chr_a > chr_b:
        chr_a, chr_b = chr_b, chr_a
        pos_a, pos_b = pos_b, pos_a
      
      bin_a = int(pos_a/sep_threshold)
      bin_b = int(pos_b/sep_threshold)
 
      key = (chr_a, chr_b, bin_a, bin_b)
      ag_data[ambig_group].append((key, pos_a, pos_b, i))
      
      if ambig_group in pos_ags: 
        pa_contacts[key].append((pos_a, pos_b, ambig_group))
      else:
        nonpa_contacts[key].append((pos_a, pos_b, ambig_group))
      
  msg = ' .. loaded {:,} contact pairs in {:,} ambiguity groups'.format(i+1, len(ag_data))
  info(msg)
 
  return ag_data, pos_ags, nonpa_contacts, pa_contacts


def get_chromo_extents(ag_data):
   
  chr_max = defaultdict(int)

  for ag in ag_data:
     for ((chr_a, chr_b, bin_a, bin_b), pos_a, pos_b, i) in ag_data[ag]:
        if pos_a > chr_max[chr_a]:
            chr_max[chr_a] = pos_a

        if pos_b > chr_max[chr_b]:
            chr_max[chr_b] = pos_b
            
            
  return dict(chr_max)


def get_unambig_contacts(ag_data, nonpa_contacts):

  unambig_contacts = {}
  unambig_pair_counts = defaultdict(int)
  chromo_bins = set()
  chromos = set()
  
  for key in nonpa_contacts:
    chr_a, chr_b, bin_a, bin_b = key
    chromos.update((chr_a, chr_b))
    vals = sorted(nonpa_contacts[key]) # Sort by 1st chromo pos
    unambig = [x[:2] for x in vals if len(ag_data[x[2]]) == 1]
    unambig_contacts[key] = unambig
    nonpa_contacts[key] = vals
             
    if unambig:
      chromo_bins.add((chr_a, bin_a))
      chromo_bins.add((chr_b, bin_b))
      unambig_pair_counts[(chr_a, chr_b)] += len(unambig)

  return sorted(chromos), chromo_bins, unambig_contacts, dict(unambig_pair_counts)



def remove_isolated_unambig(in_ncc_path, out_ncc_path, threshold=ISOLATION_THRESHOLD,
                            sep_threshold=DEFAULT_SEP_THRESHOLD, homo_trans_dens_quant=90.0,
                            keep_unambig_trans_max_chr_size=1.1e8):
  
  orig_ag_data, pa_ags, nonpa_contacts, pa_contacts = _load_bin_sort_ncc(in_ncc_path, sep_threshold)
  
  unambig_contacts = get_unambig_contacts(orig_ag_data, nonpa_contacts)[2]
  
  chr_max = get_chromo_extents(orig_ag_data)

  msg_template = ' .. Processed:{:>7,} Removed:{:>7,}'
  
  all_bins = {}
  for key in set(pa_contacts) | set(nonpa_contacts):
    all_bins[key] = pa_contacts.get(key, []) + nonpa_contacts.get(key, [])
  
  removed_ag = set()
  resolved_ag = {}
  
  msg = ' Removing inter-homologue contacts from anomalously dense regions'
  info(msg)
  
  bin_counts = {}
  for key in all_bins:
    chr_a, chr_b, bin_a, bin_b = key
    root_a = chr_a.split('.')[0]
    root_b = chr_b.split('.')[0]
    
    if root_a == root_b:
      n_contacts = len(all_bins[key])
      
      if n_contacts:
        bin_counts[key] = n_contacts
  
  if not bin_counts:
    return 0
  
  counts = list(bin_counts.values())
  upper_dens_thresh = np.percentile(counts, homo_trans_dens_quant)
  
  ag_data = dict(orig_ag_data)
  
  for ag in ag_data:
    pairs = ag_data[ag]
    n_pairs = len(pairs)
    keep = []
    
    for j, (key, pos_a, pos_b, line_idx) in enumerate(pairs):
      chr_a, chr_b, bin_a, bin_b = key
 
      if (chr_a != chr_b) and (key in bin_counts): # Homolog trans only
        count = bin_counts[key]
        
        if count < upper_dens_thresh:
          for a in range(bin_a-1, bin_a+2):
            for b in range(bin_b-1, bin_b+2):
              if a == bin_a and b == bin_b:
                continue
 
              key2 =  (chr_a, chr_b, a, b)
 
              for pos1, pos2, ag1 in all_bins.get(key2, []):
                if (abs(pos1-pos_a) < sep_threshold) and (abs(pos2-pos_b) < sep_threshold):
                  count += 1
 
                  if count > upper_dens_thresh:
                    break
              
              else:
                continue
              break
            
            else:
              continue
            break  
              
        if (count < upper_dens_thresh) and (n_pairs > 1): # ONly keep ambiguous holologous trans
          keep.append(j)
        
        elif (n_pairs == 1) and (key in unambig_contacts):
          
          contacts = []
          for pos1, pos2 in unambig_contacts[key]:
            if (pos1 != pos_a) or (pos2 != pos_b):
              contacts.append((pos1, pos1))
          
          unambig_contacts[key] = contacts
          
      else:
        keep.append(j)
    
    n_keep = len(keep)
        
    if n_keep:
      if n_keep < n_pairs:
        new_pairs = [pairs[j] for j in keep]
        ag_data[ag] = new_pairs
      
        if n_keep > 1: # Still ambiguous
          resolved_ag[ag] = [x[3] for x in new_pairs] # Line indices        
    
    else:      
      removed_ag.add(ag)
    
  msg = ' .. removed {:,}'.format(len(removed_ag))
  info(msg)
 
  msg = ' .. partly resolved {:,}'.format(len(resolved_ag))
  info(msg)

  scores = []
  it = 0

  unambig_pairs = [(g, ag_data[g][0]) for g in ag_data if len(ag_data[g]) == 1 and (g not in removed_ag)] #  and (g not in resolved_ag)]
 
  msg = ' .. keep {:,} unambiguous pairs'.format(len(unambig_pairs))
  info(msg)
  #info(' .. processing')
   
  for ag, (key, pos_a, pos_b, line_idx) in unambig_pairs:

    it += 1
    if it % 10000 == 0:
      msg = msg_template.format(it, len(removed_ag))
      info(msg, prefix='\rINFO')
 
    group_score = 0.0 
    chr_a, chr_b, bin_a, bin_b = key
      
    contacts = unambig_contacts.get(key)
 
    if contacts:
      score = _get_network_score(sep_threshold, chr_a, chr_b, pos_a, pos_b, bin_a, bin_b, unambig_contacts) # , secondary_limit=threshold) 
    else:
      score = 0.0
     
    if score < threshold:
      root_a = chr_a.split('.')[0]
      root_b = chr_b.split('.')[0]
      
      involves_small = (chr_max[chr_a] <= keep_unambig_trans_max_chr_size) or (chr_max[chr_b] <= keep_unambig_trans_max_chr_size)      
      non_homolog_trans = root_a != root_b
      
      if not (non_homolog_trans and involves_small): # Do not remove unambig from normal trans when to small chromo
        removed_ag.add(ag)
  
  _write_ambig_filtered_ncc(in_ncc_path, out_ncc_path, orig_ag_data,
                            resolved_ag=resolved_ag, removed_ag=removed_ag)

  msg = msg_template.format(it, len(removed_ag))     
  info(msg)
  
  n_non_isol = len(ag_data)
  info(' .. keeping {:,} non-isolated contacts'.format(n_non_isol))
  
  return n_non_isol
  
  
def network_filter_ambig(in_ncc_path, trans_close, primary_weight, sep_threshold, score_threshold,
                         remove_isolated=False, remove_pos_ambig=False, partition_bin_size=1e6,
                         partition_bin_sigma=9, trans_relay_percentile=5.0, max_pairs=16):
  
  ag_data, pos_ags, nonpa_contacts, pa_contacts = _load_bin_sort_ncc(in_ncc_path, sep_threshold)
  chromos, chromo_bins, unambig_contacts, unambig_pair_counts = get_unambig_contacts(ag_data, nonpa_contacts)
  
  missing_cis = set()
  
  for i, chr_a in enumerate(chromos):
    for chr_b in chromos[i:]:
      n_pair = unambig_pair_counts.get((chr_a, chr_b), 0)
      
      if chr_a == chr_b:
        if n_pair < primary_weight:
          missing_cis.add((chr_a, chr_b))
          msg = f'Chromosome {chr_a} is missing!'
          warn(msg)

  msg_template = ' .. processed:{:>7,} resolved:{:>7,} step_time:{:.3f} s time taken:{:5.2f} s'
  
  if trans_close: # Filter impossible, very weak unambiguous trans chr pairs
    trans_vals = sorted([x for x in trans_close.values() if x])
    min_trans_relay = np.percentile(trans_vals, trans_relay_percentile)
    excluded_homolog_unambig = get_homolog_partition_exclusions(chromos, ag_data, unambig_contacts, unambig_pair_counts,
                                                                sep_threshold, min_trans_relay,
                                                                partition_bin_size, partition_bin_sigma)
  
  else:
    trans_close = {}
    excluded_homolog_unambig = {}
    min_trans_relay = 0
  
  def _unambig_is_excluded(chr_a, chr_b, pos_a, pos_b):
  
    excluded = False
 
    if (chr_a in excluded_homolog_unambig) and (chr_b in excluded_homolog_unambig[chr_a]):
      ebin = int(pos_a/partition_bin_size)
      excluded = ebin in excluded_homolog_unambig[chr_a][chr_b]
 
    if (not excluded) and (chr_b in excluded_homolog_unambig) and (chr_a in excluded_homolog_unambig[chr_b]):
      ebin = int(pos_b/partition_bin_size)
      excluded = ebin in excluded_homolog_unambig[chr_b][chr_a]
    
    return excluded
       
  resolved_ag = {}
  removed_ag = set()
  start_time = t0 = time()
    
  for j, ag in enumerate(ag_data):
  
    if j % 10000 == 0:
      t1 = time()
      msg = msg_template.format(j, len(resolved_ag), t1-t0, t1-start_time)
      t0 = t1
      info(msg, prefix='\rINFO')
    
    if ag in removed_ag:
      continue
    
    pairs = ag_data[ag]
 
    if len(pairs) == 1:
      (chr_a, chr_b, bin_a, bin_b), pos_a, pos_b, line_idx = pairs[0]
      if _unambig_is_excluded(chr_a, chr_b, pos_a, pos_b): # Wrong side of a homologue partition
        removed_ag.add(ag)

      continue
      
    if len(pairs) > max_pairs:
      removed_ag.add(ag)
      continue

    poss_pairs = []
    for key, pos_a, pos_b, line_idx in pairs:
      chr_pair = tuple(key[:2])
      
      chr_a, chr_b = chr_pair
      
      if chr_a in missing_cis:
        continue
        
      if chr_b in missing_cis:
        continue
      
      # cis or significant trans  
      if (chr_pair not in trans_close) or (trans_close[chr_pair] >= min_trans_relay): 
        poss_pairs.append((key, pos_a, pos_b, line_idx))
    
    n_poss = len(poss_pairs)
    if n_poss < len(pairs): # A chromo pair can be excluded
      if n_poss == 0:
        removed_ag.add(ag)
    
      elif n_poss == 1: # Only one sensible trans possibility
        key, pos_a, pos_b, line_idx = poss_pairs[0]
        chr_a, chr_b, bin_a, bin_b = key

        if nonpa_contacts.get(key) and (_get_network_score(sep_threshold, chr_a, chr_b, pos_a, pos_b, bin_a, bin_b, unambig_contacts) > 0.0):
          if _unambig_is_excluded(chr_a, chr_b, pos_a, pos_b): # Wrong side of a homologue partition
            removed_ag.add(ag)
          else:
            resolved_ag[ag] = (line_idx,)
        else:
          removed_ag.add(ag)
        
        continue
      
      else:
        resolved_ag[ag] = tuple([x[3] for x in poss_pairs]) # Carry on, might be refined further

    line_indices = []
    m = max(4, len(poss_pairs))
    group_scores = [0.0] * m
    unambig_exclusions = [None] * m
    
    # Homolog partition exclusions prevent unambiguous in exclusion zone
    # - Isolated unambiguous removed
    # - Cannot become unambiguous
    
    for p, (key, pos_a, pos_b, line_idx) in enumerate(poss_pairs):
      line_indices.append(line_idx)
      chr_a, chr_b, bin_a, bin_b = key
      contacts = nonpa_contacts.get(key)
     
      if chr_a in missing_cis:
        score = 0.0

      elif chr_b in missing_cis:
        score = 0.0
 
      elif contacts:
        score = _get_network_score(sep_threshold, chr_a, chr_b, pos_a, pos_b, bin_a, bin_b, unambig_contacts)
 
      else:
        score = 0.0
      
      unambig_exclusions[p] = _unambig_is_excluded(chr_a, chr_b, pos_a, pos_b)
      group_scores[p] = score

    a, b, c, d, *e = np.argsort(group_scores)[::-1]
    
    if group_scores[a] < 1.0: # Best is isolated
      if remove_isolated or unambig_exclusions[a]:
        if ag in resolved_ag:
          del resolved_ag[ag]
          
        removed_ag.add(ag)
     
    elif group_scores[a] > score_threshold * group_scores[b]:
      if unambig_exclusions[a]:
        if group_scores[b] < 1.0: # No alternative
          if ag in resolved_ag:
            del resolved_ag[ag]
          
          removed_ag.add(ag)
          
        else:
          resolved_ag[ag] = (line_indices[a], line_indices[b])
        
      else: # Becomes unambig
        resolved_ag[ag] = (line_indices[a],)
      
    elif group_scores[b] > score_threshold * group_scores[c]: # First two were close
      resolved_ag[ag] = (line_indices[a], line_indices[b])

    elif group_scores[c] > score_threshold * group_scores[d]: # First three were close
      resolved_ag[ag] = (line_indices[a], line_indices[b], line_indices[c])
    
  if remove_pos_ambig: 
    removed_ag.update(pos_ags)
    
  return ag_data, resolved_ag, removed_ag


def get_homolog_partition_exclusions(chromos, ag_data, unambig_contacts, unambig_pair_counts,
                                     sep_threshold, min_homolog_trans_count,
                                     partition_bin_size, partition_bin_sigma):

  excluded_unambig = {}
  chr_roots = sorted(set([x.split('.')[0] for x in chromos]))
  homolog_pairs = [(f'{x}.a', f'{x}.b') for x in chr_roots]
  homolog_pairs = [(a,b) for a, b in homolog_pairs if ((a in chromos) and (b in chromos))] # A homologue can be missing, esp. for chrX  
  chr_max = get_chromo_extents(ag_data)  
  width = 1 + 8 * partition_bin_sigma
  kernel = signal.windows.gaussian(width, std=partition_bin_sigma)
  
  from matplotlib import pyplot as plt
  
  for chr_1 in chromos:
      pos_max_1 = chr_max[chr_1]
      bin_max_1 = int(ceil(pos_max_1/sep_threshold))
      nbins_1 = 1 + int(ceil(pos_max_1/partition_bin_size))
      excluded_unambig[chr_1] = defaultdict(set)
  
      for chr_2a, chr_2b in homolog_pairs:

          pos_max_2a = chr_max[chr_2a]
          pos_max_2b = chr_max[chr_2b]
          bin_max_2a = int(ceil(pos_max_2a/sep_threshold))
          bin_max_2b = int(ceil(pos_max_2b/sep_threshold))
  
          pair_key_a = tuple(sorted([chr_1, chr_2a]))
          pair_key_b = tuple(sorted([chr_1, chr_2b]))
   
          counts_a = unambig_pair_counts.get(pair_key_a, 0)
          counts_b = unambig_pair_counts.get(pair_key_b, 0)
  
          if (counts_a > min_homolog_trans_count) and (counts_b > min_homolog_trans_count):
              influence = np.zeros(nbins_1)
              
              # set to +1 for chr_2a unambig contacts
              # set to -1 for chr_2b unambig contacts
              
              info(f' .. partitioning {chr_1} vs {chr_2a}[{counts_a}]:{chr_2b}[{counts_b}]')
              
              for bin1 in range(bin_max_1):
                  for bin2a in range(bin_max_2a):
                      key = (chr_1, chr_2a, bin1, bin2a)
                      for pos1, pos2a in unambig_contacts.get(key, []):
                          ebin = int(pos1/partition_bin_size) # This bin on chr_1
                          influence[ebin] += 1.0

                      key = (chr_2a, chr_1, bin2a, bin1)
                      for pos2a, pos1 in unambig_contacts.get(key, []):
                          ebin = int(pos1/partition_bin_size)
                          influence[ebin] += 1.0

                  for bin2b in range(bin_max_2b):
                      key = (chr_1, chr_2b, bin1, bin2b)
                      for pos1, pos2b in unambig_contacts.get(key, []):
                          ebin = int(pos1/partition_bin_size)
                          influence[ebin] -= 1.0

                      key = (chr_2b, chr_1, bin2b, bin1)
                      for pos2b, pos1 in unambig_contacts.get(key, []):
                          ebin = int(pos1/partition_bin_size)
                          influence[ebin] -= 1.0
              
              # # Plots to debug & visualise homologue exclusion partitioning
              #plt.plot(influence)
              influence = np.convolve(influence, kernel, 'full')
              p = len(kernel)//2
              influence = influence[p:p+nbins_1]
              allowed_bins_a = list(np.flatnonzero(influence > 0.0))
              allowed_bins_b = list(np.flatnonzero(influence < 0.0))
              
              #plt.plot(influence)
              #plt.show()
              
              # Mutual exclusion between homolog chromos to other chr
              excluded_unambig[chr_1][chr_2a].update(allowed_bins_b)
              excluded_unambig[chr_1][chr_2b].update(allowed_bins_a)
  
  return excluded_unambig
                           

def get_trans_closeness(ncc_path, ag_data, sep_threshold, primary_weight):

  nonpa_contacts = _load_bin_sort_ncc(ncc_path, sep_threshold)[2]
  chromos, chromo_bins, unambig_contacts, unambig_pair_counts = get_unambig_contacts(ag_data, nonpa_contacts)    
 
  trans_close = {}  
  for i, chr_a in enumerate(chromos):
    for chr_b in chromos[i:]:

      if chr_a != chr_b:
        n_pair = primary_weight * unambig_pair_counts.get((chr_a, chr_b), 0)
        
        for chr_c in chromos:
          if chr_c in (chr_a, chr_b):
            continue
 
          a = unambig_pair_counts.get((chr_a, chr_c), 0)
          b = unambig_pair_counts.get((chr_b, chr_c), 0)
          c = unambig_pair_counts.get((chr_c, chr_a), 0)
          d = unambig_pair_counts.get((chr_c, chr_b), 0)
 
          n_pair += np.sqrt((a + c) * (b + d))
 
        trans_close[(chr_a, chr_b)] = n_pair
        trans_close[(chr_b, chr_a)] = n_pair
  
  return trans_close
  
  
def resolve_contacts(in_ncc_path, out_ncc_path, temp_ncc_path1, sep_threshold, remove_isolated=True,
                     score_threshold=2.0, remove_pos_ambig=False, primary_weight=5,
                     trans_relay_percentile=5.0, partition_bin_size=1e6, partition_bin_sigma=9):
  
  pass_args = (primary_weight, sep_threshold, score_threshold, remove_isolated,
               remove_pos_ambig, partition_bin_size, partition_bin_sigma, trans_relay_percentile)
               
  msg_template = ' .. processed:{:>7,} resolved:{:>7,}'
  
  # Do a light disambiguation initially, which gives better trans relay stats
  # - min_trans_relay comes from stats

  info(' Primary ambiguity filtering with network scores')
  
  ag_data, resolved_ag, removed_ag = network_filter_ambig(in_ncc_path, None, *pass_args[:3]) 
  info(msg_template.format(len(ag_data), len(resolved_ag)))  

  info(' Filtering weak trans chromsosome pairs')
  
  _write_ambig_filtered_ncc(in_ncc_path, temp_ncc_path1, ag_data, resolved_ag, removed_ag)  
  trans_close = get_trans_closeness(temp_ncc_path1, ag_data, sep_threshold, primary_weight)
  
  info(' Secondary ambiguity filtering with network scores')

  ag_data, resolved_ag, removed_ag = network_filter_ambig(in_ncc_path, trans_close, *pass_args)  
  
  _write_ambig_filtered_ncc(in_ncc_path, out_ncc_path, ag_data, resolved_ag, removed_ag)  
  info(msg_template.format(len(ag_data), len(resolved_ag)))  
 

def sc_hic_disambiguate(ncc_file_paths, sep_threshold=DEFAULT_SEP_THRESHOLD, score_threshold=DEFAULT_SCORE_THRESHOLD,
                        min_contacts=DEFAULT_MIN_CONTACTS, out_suffix=DEFAULT_OUT_SUFFIX,
                        keep_unambig_trans_max_chr_size=DEFAULT_SMALL_CHR_SIZE_LIMIT,
                        keep_intermed=False,  intermed_maps=False):

  from tools.contact_map import contact_map
  from hic_core.nuc_process import SESSION_KEY
  
  info('Resolving ambiguous single-cell Hi-C contacts')
    
  for ncc_file_path in ncc_file_paths:
    
    file_root = os.path.splitext(ncc_file_path)[0] + out_suffix
    
    clean_ncc_path = file_root + '_' + SESSION_KEY + '_stage1_clean.ncc'
    temp_ncc_path1 = file_root + '_' + SESSION_KEY + '_stage2_filter.ncc'
    out_ncc_path = file_root + '.ncc'
    
    n = remove_isolated_unambig(ncc_file_path, clean_ncc_path,
                                sep_threshold=sep_threshold,
                                keep_unambig_trans_max_chr_size=keep_unambig_trans_max_chr_size)
    
    if n < min_contacts:
      msg = 'File {} will not be processed further. Number of non-isolated contacts: {:,} is smaller than minimum required: {:,}.'
      warn(msg.format(ncc_file_path, n, min_contacts))

      if intermed_maps:
        clean_pdf = file_root + '_stage1_clean_map.pdf'
        contact_map([clean_ncc_path], clean_pdf, bin_size=5e3, show_chromos=None,
                     no_separate_cis=True, is_single_cell=True)
      
      if keep_intermed:
        shutil.move(clean_ncc_path, file_root + '_stage1_clean.ncc')
      else:
        os.unlink(clean_ncc_path)      
      
    else:
      if intermed_maps:
        clean_pdf = file_root + '_stage1_clean_map.pdf'
        contact_map([clean_ncc_path], clean_pdf, bin_size=5e3, show_chromos=None,
                     no_separate_cis=True, is_single_cell=True)
        
      resolve_contacts(clean_ncc_path, out_ncc_path, temp_ncc_path1, sep_threshold, remove_isolated=True, score_threshold=score_threshold)

      if intermed_maps:
        filter_pdf = file_root + '_stage2_filter_map.pdf'
        contact_map([temp_ncc_path1], filter_pdf, bin_size=5e3, show_chromos=None,
                     no_separate_cis=True, is_single_cell=True)

      if keep_intermed:
        shutil.move(clean_ncc_path, file_root + '_stage1_clean.ncc')
        shutil.move(temp_ncc_path1, file_root + '_stage2_filter.ncc')
      else:
        os.unlink(clean_ncc_path)
        os.unlink(temp_ncc_path1)
      
      out_main_pdf = file_root + '_map.pdf'
      contact_map([out_ncc_path], out_main_pdf, bin_size=5e3, show_chromos=None,
                   no_separate_cis=True, is_single_cell=True)

  
  
def main(argv=None):
  from argparse import ArgumentParser

  if argv is None:
    argv = sys.argv[1:]

  epilog = 'For further help email tstevens@mrclmb.ac.uk or wb104@cam.ac.uk'

  arg_parse = ArgumentParser(prog=PROG_NAME, description=DESCRIPTION,
                             epilog=epilog, prefix_chars='-', add_help=True)

  arg_parse.add_argument(nargs='+', metavar='NCC_FILE', dest='i',
                         help='One or more single-cell Hi-C contact files in NCC format. Wildcards accepted.')

  arg_parse.add_argument('-s', '--seq-sep-max', metavar='SEQ_SEP', default=DEFAULT_SEP_THRESHOLD/1e6, type=int, dest='s',
                         help='Sequence separation threshold (in Mb) within which contacts may be considered as supporting. ' \
                              'Default: {:,} Mb'.format(DEFAULT_SEP_THRESHOLD/1e6))

  arg_parse.add_argument('-t', '--score-threshold', metavar='MIN_SCORE', default=DEFAULT_SCORE_THRESHOLD, type=float, dest='t',
                         help='Network support score threshold above which unambiguous contacts help resolve ambiguous ones. ' \
                              'A higher value means more strict. Default: {:,} (cautious)'.format(DEFAULT_SCORE_THRESHOLD))

  arg_parse.add_argument('-n', '--min-num-contacts', metavar='MIN_CONTACTS', default=DEFAULT_MIN_CONTACTS, type=int, dest='n',
                         help='Minimum number of non-isolated contacts required overall for disambiguation. Default: {:,}'.format(DEFAULT_MIN_CONTACTS))

  arg_parse.add_argument('-o', '--out-suffix', metavar='OUT_SUFFIX', default=DEFAULT_OUT_SUFFIX, dest='o',
                         help='File suffix to add to input NCC file name(s) to make the filtered output NCC file name and PDF contact map file name,' \
                              ' e.g from "Cell_1.ncc" to "Cell_1_disambig.ncc" and "Cell_1_disambig_map.pdf". Default: {}'.format(DEFAULT_OUT_SUFFIX))

  arg_parse.add_argument('-z', '--small-chr_size', metavar='SEQ_SIZE', default=DEFAULT_SMALL_CHR_SIZE_LIMIT/1e6, type=int, dest='z',
                         help='Only chromosomes larger than this size can have isolated unambiguous (non-homologue) trans contacts removed. ' \
                              'Default: {:,} Mb'.format(DEFAULT_SMALL_CHR_SIZE_LIMIT/1e6))

  arg_parse.add_argument('-k', '--keep-intermed', default=False, action='store_true', dest='k',
                         help='Keep NCC contact files from intermediate processing stages; else only final contacts will be kept.')

  arg_parse.add_argument('-m', '--intermed-maps', default=False, action='store_true', dest='m',
                         help='Make contact maps for intermediate disambiguation stages; else only a final contact map is made.')

  args = vars(arg_parse.parse_args(argv))

  ncc_paths = args['i']
  sep_threshold = args['s']
  score_threshold = args['t']
  min_contacts = args['n']
  out_suffix = args['o']
  intermed_maps = args['m']
  keep_intermed = args['k']
  keep_unambig_trans_max_chr_size = args['z']
  
  if not ncc_paths:
    msg = 'No NCC format contact files specified'
    fatal(msg)

  for ncc_path in ncc_paths:
    if not os.path.exists(ncc_path):
      msg = 'NCC file {} does not exist'.format(ncc_path)
      fatal(msg)
  
  if not (1.0 <= sep_threshold <= 100.0):
    msg = 'Sequence separation threshold must be between 1 and 100 Mb'
    fatal(msg)
  
  if not (1.0 <= score_threshold <= 100.0):
    msg = 'Score threshold must be between 1.0 and 100.0'
    fatal(msg)
  
  sep_threshold *= 1e6
  keep_unambig_trans_max_chr_size *= 1e6
  
  sc_hic_disambiguate(ncc_paths, sep_threshold, score_threshold, min_contacts, out_suffix,
                      keep_unambig_trans_max_chr_size, keep_intermed, intermed_maps)


if __name__ == '__main__':
  main_dir = os.path.dirname(os.path.dirname(__file__))
  sys.path.insert(0, main_dir)
  sys.path.insert(0, os.path.join(main_dir, 'nuc_tools')) 
  main()


