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
import os, sys, shutil
from math import log, exp
from time import time
from collections import defaultdict

import numpy as np

from hic_core.nuc_process import info, warn, fatal

PROG_NAME = 'single_cell_hic_disambiguate'
VERSION = '1.0.0'
DESCRIPTION = 'An extension toool for nuc_process that attempts to resolve ambiguous (multi-mapping) contacts in single-cell Hi-C data. Operates on NCC format files.'

MIN_SCORE_SEQ_SEP = 100
ISOLATION_THRESHOLD = 0.01

DEFAULT_SEP_THRESHOLD  = 4e7
DEFAULT_SCORE_THRESHOLD = 10.0
DEFAULT_MIN_CONTACTS = int(1e5)
DEFAULT_OUT_SUFFIX = '_disambig'

def _get_network_score(sep_threshold, chr_a, chr_b, pos_a, pos_b, bin_a, bin_b,
                       unambig_bins, chromo_bins, primary_limit=8.0, min_sep=MIN_SCORE_SEQ_SEP):
  """
  primary_limit is an early-stopping heristic to improve performance; should not be changed
  """
  sep_scale = 0.5 * sep_threshold
  plim = 10.0 ** primary_limit
  score_l = 1.0
  score_u = 1.0
  unambig_bins_get = unambig_bins.get
  
  a_bins = (bin_a-1, bin_a, bin_a+1)
  b_bins = (bin_b-1, bin_b, bin_b+1)
  
  for a in a_bins:
    rev = a < bin_a
    peri_a = a != bin_a
    
    for b in b_bins:
      key2 = (chr_a, chr_b, a, b)
      
      unambig = unambig_bins_get(key2)
      
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
        if int(float(ambig_code)) > 0:
          ambig_group += 1
      else:
        ambig_group = int(ambig_code)    
      
      sz = len(ag_data[ambig_group])
      
      if ambig_group in resolved_ag:    
        keep = 1 if i in resolved_ag[ambig_group] else 0 # Else inactivate
        
        if ambig_group in seen_res_group:
          row[12] = '0.%d' % (keep,)
        else:
          seen_res_group.add(ambig_group)
          row[12] = '%d.%d'  % (sz,keep)
                    
        line = ' '.join(row)+'\n'

      elif ambig_group in removed_ag: # Inactivate
        if ambig_group in seen_res_group:
          row[12] = '0.0'
        else:
          seen_res_group.add(ambig_group)
          row[12] = '%d.0'  % (sz,)
           
        line = ' '.join(row)+'\n'
        
      write(line)    
      
      n += 1

  msg = " .. written {:,} of {:,} lines to {}".format(n, i+1, out_ncc_path)
  
  info(msg)
  
      
def _load_bin_sort_ncc(in_ncc_path, sep_threshold):
  
  msg = ' .. reading contacts from %s' % in_ncc_path
  info(msg)

  ag_data = defaultdict(list)
  nonambig_bins = defaultdict(list)
  ag_positions = defaultdict(set)
  pos_ambig = set()
  
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
  
  n_pos_ag = 0
  for ag in ag_positions:
    chr_pos = ag_positions[ag]
    
    if len(chr_pos) > 4:
      pos_ambig.add(ag)
      
    else:
      chr_roots = set([x[0].split('.')[0] for x in chr_pos])
      
      if len(chr_roots) > 2:
        pos_ambig.add(ag)
  
  msg = ' .. found {:,} positional ambiguity groups from {:,} total'.format(len(pos_ambig), ambig_group)
  info(msg)
    
  # Bin and group
  ambig_bins = defaultdict(list)
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
      
      
      if ambig_group in pos_ambig: 
        ambig_bins[key].append((pos_a, pos_b, ambig_group))
      else:
        nonambig_bins[key].append((pos_a, pos_b, ambig_group))
      
  msg = ' .. loaded {:,} contact pairs in {:,} ambiguity groups'.format(i+1, len(ag_data))
  info(msg)
 
  msg = ' .. sorting data'
  info(msg)

  unambig_bins = {}
  chromo_pair_counts = defaultdict(int)
  chromo_bins = set()
  chromos = set()
  
  for key in nonambig_bins:
    chr_a, chr_b, bin_a, bin_b = key
    chromos.update((chr_a, chr_b))
    vals = sorted(nonambig_bins[key]) # Sort by 1st chromo pos
    unambig = [x[:2] for x in vals if len(ag_data[x[2]]) == 1]
 
    unambig_bins[key] = unambig
    nonambig_bins[key] = vals
             
    if unambig:
      chromo_bins.add((chr_a, bin_a))
      chromo_bins.add((chr_b, bin_b))
      
      chromo_pair_counts[(chr_a, chr_b)] += len(unambig)

  return sorted(chromos), ag_data, pos_ambig, chromo_bins, nonambig_bins, unambig_bins, ambig_bins, dict(chromo_pair_counts)


def remove_isolated_unambig(in_ncc_path, out_ncc_path, threshold=ISOLATION_THRESHOLD, sep_threshold=DEFAULT_SEP_THRESHOLD, homo_trans_dens_quant=90.0):
  
  chromos, ag_data, pos_ambig, chromo_bins, nonambig_bins, unambig_bins, ambig_bins, chromo_pair_counts = _load_bin_sort_ncc(in_ncc_path, sep_threshold)
 
  msg_template = ' .. processed:{:>7,} removed:{:>7,}'
  
  all_bins = {}
  for key in set(ambig_bins) | set(nonambig_bins):
    all_bins[key] = ambig_bins.get(key, []) + nonambig_bins.get(key, [])
  
  removed_ag = set()
  resolved_ag = {}
  
  msg = ' Removing inter-homologue contacts from anomolously dense regions'
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
        
        elif (n_pairs == 1) and (key in unambig_bins):
          
          contacts = []
          for pos1, pos2 in unambig_bins[key]:
            if (pos1 != pos_a) or (pos2 != pos_b):
              contacts.append((pos1, pos1))
          
          unambig_bins[key] = contacts
          
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
      
    contacts = unambig_bins.get(key)
 
    if contacts:
      score = _get_network_score(sep_threshold, chr_a, chr_b, pos_a, pos_b, bin_a, bin_b, unambig_bins, chromo_bins) # , secondary_limit=threshold) 
    else:
      score = 0.0
     
    if score < threshold:
      removed_ag.add(ag)
  
  
  _write_ambig_filtered_ncc(in_ncc_path, out_ncc_path, ag_data, resolved_ag=resolved_ag, removed_ag=removed_ag)

  msg = msg_template.format(it, len(removed_ag))     
  info(msg)
  
  n_non_isol = len(ag_data)
  info(' .. keeping {:,} non-isolated contacts'.format(n_non_isol))
  
  return n_non_isol
  
  
def network_filter_ambig(ag_data, missing_cis, trans_close, nonambig_bins, unambig_bins, chromo_bins,
                         sep_threshold, score_threshold, min_trans_relay=5, remove_isolated=False, max_pairs=16):

  msg_template = ' .. processed:{:>7,} resolved:{:>7,} step_time:{:.3f} s time taken:{:5.2f} s'
  
  resolved_ag = {}
  removed_ag = set()
  
  scores = [[], [], [], []]
  start_time = t0 = time()
  
  
  qc = []
  #info(' .. processing')
    
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
      continue
      
    if len(pairs) > max_pairs:
      removed_ag.add(ag)
      continue

    poss_pairs = []
    for key, pos_a, pos_b, line_idx in pairs:
      chr_pair = tuple(key[:2])
      
      if chr_pair in missing_cis:
        continue
        
      if (chr_pair not in trans_close) or (trans_close[chr_pair] >= min_trans_relay): 
        poss_pairs.append((key, pos_a, pos_b, line_idx))
    
    n_poss = len(poss_pairs)
    if poss_pairs and n_poss < len(pairs): # A chromo pair can be excluded
      if len(poss_pairs) == 1: # Only one sensible trans possibility
         key, pos_a, pos_b, line_idx = poss_pairs[0]
         chr_a, chr_b, bin_a, bin_b = key
          
         if nonambig_bins.get(key) and _get_network_score(sep_threshold, chr_a, chr_b, pos_a, pos_b, bin_a, bin_b, unambig_bins, chromo_bins):
           resolved_ag[ag] = (line_idx,)
         else:
           removed_ag.add(ag)
        
         continue
      
      else:
        resolved_ag[ag] = tuple([x[3] for x in poss_pairs]) # Might be refined further

    line_indices = []
    group_scores = [0.0] * max(4, len(pairs))
    
    for p, (key, pos_a, pos_b, line_idx) in enumerate(pairs):
      line_indices.append(line_idx)
      chr_a, chr_b, bin_a, bin_b = key
      contacts = nonambig_bins.get(key)
     
      if (chr_a, chr_b) in missing_cis:
        score = 0.0
 
      elif contacts:
        score = _get_network_score(sep_threshold, chr_a, chr_b, pos_a, pos_b, bin_a, bin_b, unambig_bins, chromo_bins)
 
        if score and chr_a != chr_b:
          qc.append(score)
 
      else:
        score = 0.0
        
      group_scores[p] = score

    a, b, c, d, *e = np.argsort(group_scores)[::-1]
    
    if group_scores[a] < 1.0: # Best is isolated
      if remove_isolated:
        if ag in resolved_ag:
          del resolved_ag[ag]
          
        removed_ag.add(ag)
     
    elif group_scores[a]  >  score_threshold * group_scores[b]:
      resolved_ag[ag] = (line_indices[a],)
      
    elif group_scores[b]  >  score_threshold * group_scores[c]: # First two were close
      resolved_ag[ag] = (line_indices[a], line_indices[b])

    elif group_scores[c]  >  score_threshold * group_scores[d]: # First three were close
      resolved_ag[ag] = (line_indices[a], line_indices[b], line_indices[c])
    
  return resolved_ag, removed_ag

  
def resolve_contacts(in_ncc_path, out_ncc_path, temp_ncc_path, sep_threshold, remove_isolated=True, score_threshold=2.0,
                     remove_pos_ambig=False, primary_weight=5, trans_relay_percentile=5.0):
  
  msg_template = ' .. processed:{:>7,} resolved:{:>7,}'

  msg = ' Reading contact data'
  info(msg)
  
  chromos, ag_data, pos_ambig, chromo_bins, nonambig_bins, unambig_bins, ambig_bins, chromo_pair_counts = _load_bin_sort_ncc(in_ncc_path, sep_threshold)

  missing_cis = set()
  
  for i, chr_a in enumerate(chromos):
    for chr_b in chromos[i:]:
      n_pair = chromo_pair_counts.get((chr_a, chr_b), 0)
      
      if chr_a == chr_b:
        if n_pair < primary_weight:
          missing_cis.add((chr_a, chr_b))
          msg = f'Chromosome {chr_a} is missing!'
          warn(msg)
  
  # Do a light disambiguation initially, which gives better trans relay stats
  # - min_trans_relay comes from stats

  msg = ' Primary ambiguity filtering with network scores'
  info(msg)
  
  resolved_ag, removed_ag = network_filter_ambig(ag_data, missing_cis, {}, nonambig_bins, unambig_bins,
                                                 chromo_bins, sep_threshold, score_threshold)

  _write_ambig_filtered_ncc(in_ncc_path, temp_ncc_path, ag_data, resolved_ag, removed_ag)
  
  chromos, ag_data, pos_ambig, chromo_bins, nonambig_bins, unambig_bins, ambig_bins, chromo_pair_counts = _load_bin_sort_ncc(temp_ncc_path, sep_threshold)
  
  info(msg_template.format(len(ag_data), len(resolved_ag)))  
  
  trans_close = {}  
  for i, chr_a in enumerate(chromos):
    for chr_b in chromos[i:]:

      if chr_a != chr_b:
        n_pair = primary_weight * chromo_pair_counts.get((chr_a, chr_b), 0)
        
        for chr_c in chromos:
          if chr_c in (chr_a, chr_b):
            continue
 
          a = chromo_pair_counts.get((chr_a, chr_c), 0)
          b = chromo_pair_counts.get((chr_b, chr_c), 0)
          c = chromo_pair_counts.get((chr_c, chr_a), 0)
          d = chromo_pair_counts.get((chr_c, chr_b), 0)
 
          n_pair += np.sqrt((a + c) * (b + d))
 
        trans_close[(chr_a, chr_b)] = n_pair
        trans_close[(chr_b, chr_a)] = n_pair

  
  trans_vals = sorted([x for x in trans_close.values() if x])

  min_trans_relay = np.percentile(trans_vals, trans_relay_percentile)

  msg = ' Secondary ambiguity filtering with network scores'
  info(msg)

  chromos, ag_data, pos_ambig, chromo_bins, nonambig_bins, unambig_bins, ambig_bins, chromo_pair_counts = _load_bin_sort_ncc(in_ncc_path, sep_threshold)

  resolved_ag, removed_ag = network_filter_ambig(ag_data, missing_cis, trans_close, nonambig_bins, unambig_bins,
                                                 chromo_bins, sep_threshold, score_threshold, min_trans_relay, remove_isolated)  

  if remove_pos_ambig: 
    removed_ag.update(pos_ambig)
  
  _write_ambig_filtered_ncc(in_ncc_path, out_ncc_path, ag_data, resolved_ag, removed_ag)
  
  info(msg_template.format(len(ag_data), len(resolved_ag)))  


def sc_hic_disambiguate(ncc_file_paths, sep_threshold=DEFAULT_SEP_THRESHOLD, score_threshold=DEFAULT_SCORE_THRESHOLD,
                        min_contacts=DEFAULT_MIN_CONTACTS, out_suffix=DEFAULT_OUT_SUFFIX,
                        keep_intermed=False,  intermed_maps=False):

  from tools.contact_map import contact_map
  
  from hic_core.nuc_process import SESSION_KEY
  
  info('Resolving ambiguous single-cell Hi-C contacts')
    
  for ncc_file_path in ncc_file_paths:
    
    file_root = os.path.splitext(ncc_file_path)[0] + out_suffix
    
    clean_ncc_path = file_root + '_' + SESSION_KEY + '_stage1_clean.ncc'
    temp_ncc_path = file_root + '_' + SESSION_KEY + '_stage2_filter.ncc'
    out_ncc_path = file_root + '.ncc'
    
    n = remove_isolated_unambig(ncc_file_path, clean_ncc_path, sep_threshold=sep_threshold)
    
    if n < min_contacts:
      msg = 'File {} will not be processed further. Number of non-isolated contacts: {:,} is smaller than minimum required: {:,}.'
      warn(msg.format(ncc_file_path, n, min_contacts))

      if intermed_maps:
        clean_pdf = file_root + '_stage1_clean_map.pdf'
        contact_map([clean_ncc_path], clean_pdf, bin_size=None, show_chromos=None,
                     no_separate_cis=True, is_single_cell=True)
      
      if keep_intermed:
        shutil.move(clean_ncc_path, file_root + '_stage1_clean.ncc')
      else:
        os.unlink(clean_ncc_path)      
      
    else:
      if intermed_maps:
        clean_pdf = file_root + '_stage1_clean_map.pdf'
        contact_map([clean_ncc_path], clean_pdf, bin_size=None, show_chromos=None,
                     no_separate_cis=True, is_single_cell=True)
        
      resolve_contacts(clean_ncc_path, out_ncc_path, temp_ncc_path, sep_threshold, remove_isolated=True, score_threshold=score_threshold)

      if intermed_maps:
        filter_pdf = file_root + '_stage2_filter_map.pdf'
        contact_map([temp_ncc_path], filter_pdf, bin_size=None, show_chromos=None,
                     no_separate_cis=True, is_single_cell=True)
      
      if keep_intermed:
        shutil.move(clean_ncc_path, file_root + '_stage1_clean.ncc')
        shutil.move(temp_ncc_path, file_root + '_stage2_filter.ncc')
      else:
        os.unlink(clean_ncc_path)
        os.unlink(temp_ncc_path)
      
      out_main_pdf = file_root + '_map.pdf'
      contact_map([out_ncc_path], out_main_pdf, bin_size=None, show_chromos=None,
                   no_separate_cis=True, is_single_cell=True)

  
  
def main(argv=None):
  from argparse import ArgumentParser

  if argv is None:
    argv = sys.argv[1:]

  epilog = 'For further help email tjs23@cam.ac.uk or wb104@cam.ac.uk'
  arg_parse = ArgumentParser(prog=PROG_NAME, description=DESCRIPTION,
                             epilog=epilog, prefix_chars='-', add_help=True)

  arg_parse.add_argument(nargs='+', metavar='FASTQ_FILE', dest='i',
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
  
  sc_hic_disambiguate(ncc_paths, sep_threshold, score_threshold, min_contacts, out_suffix, keep_intermed, intermed_maps)

if __name__ == '__main__':
  main_dir = os.path.dirname(os.path.dirname(__file__))
  sys.path.insert(0, main_dir)
  sys.path.insert(0, os.path.join(main_dir, 'nuc_tools')) 
  main()


