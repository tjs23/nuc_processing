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

PROG_NAME = 'nuc_sequence_names'
VERSION = '1.0.0'
DESCRIPTION = 'Fetch human readable chromosome/sequence names from NCBI genome build names'
DEFAULT_MIN_SIZE = 100

import sys
from collections import defaultdict

try:
  from urllib2 import HTTPError

except ImportError:
  from urllib.error import HTTPError

from NucProcess import info, warn, fatal, check_regular_file

try:
  from Bio import Entrez

except ImportError:
  msg = 'BioPython module not installed or not accessible'
  fatal(msg)

def get_chromosome_names(seq_ids, email):

  name_dict = {}
  unloc_counts = defaultdict(int)
  
  Entrez.email = email
  
  for seq_id in seq_ids:
    
    if seq_id.startswith('scaffold'):
      name_dict[seq_id] = seq_id
      continue      
    
    try:
      stream_obj = Entrez.efetch(db="nuccore", id=seq_id, format='xml')
      
    except HTTPError:
      warn('NCBI Entrez seq_id lookup failed for "%s"' % seq_id)
      name_dict[seq_id] = seq_id
      continue
    
    results = Entrez.read(stream_obj)[0]
    db_id = results['GBSeq_accession-version']
    
    chr_name = None
    map_type = None
    organelle = None

    for sub_dict in results['GBSeq_feature-table'][0]['GBFeature_quals']:
      if sub_dict['GBQualifier_name'] == 'map':
        map_type = sub_dict['GBQualifier_value']
      
      if sub_dict['GBQualifier_name'] == 'chromosome':
        chr_name = sub_dict['GBQualifier_value']

      if sub_dict['GBQualifier_name'] == 'organelle':
        organelle = sub_dict['GBQualifier_value']
    
    if organelle:
      if organelle == 'mitochondrion':
        chr_name = 'chrMt'
        
      elif organelle == 'plastid':
        chr_name = 'chrPl'
        
      else:
        chr_name = 'organelle'
      
    elif map_type == 'unlocalized':
      if chr_name:
        unloc_counts[chr_name] += 1
        chr_name = 'chr%s_unloc_%d' % (chr_name, unloc_counts[chr_name])
      
      else:
        chr_name = db_id

      info('Using name "%s" for seq_id "%s"' % (chr_name, seq_id))
     
    elif chr_name:
    
      if chr_name.lower() == 'unknown':
        info('Chromosome unknown for seq_id "%s"' % seq_id)
        chr_name = db_id
    
      else:
        chr_name = 'chr' + chr_name
        info('Using name "%s" for seq_id "%s"' % (chr_name, seq_id))
      
    else:
      warn('No chromosome name found for seq_id "%s"' % seq_id)
      chr_name = db_id

    
    name_dict[seq_id] = chr_name  
  
  return name_dict  


def write_chromosome_name_file(out_path, email, fasta_paths, min_kb_size=DEFAULT_MIN_SIZE):
  
  from NucProcess import open_file_r
  min_seq_len = min_kb_size * 1000
  seq_ids = {}
  count = 0
  
  for fasta_path in fasta_paths:    
    check_regular_file(fasta_path, critical=True)
    seq_id = None
    prev_id = None
    
    with open_file_r(fasta_path) as file_obj:
      
      for line in file_obj:
        if line[0] == '>':
          seq_id = line[1:].split()[0]
          
          if seq_id.startswith('gi|'):
            seq_id = seq_id.split('|')[3]
          
          if seq_id in seq_ids:
            fatal('Sequence ID %s repeated in FASTA file %s' % (seq_id, fasta_path))
          
          if prev_id and count > min_seq_len:
            info('Found sequence ID {} length {:,}'.format(prev_id, count))
            seq_ids[prev_id] = count
            
          count = 0
          prev_id = seq_id
        
        else:
          count += len(line) -1
  
    if prev_id and count > min_seq_len:
      info('Found sequence ID {} length {:,}'.format(prev_id, count))
      seq_ids[prev_id] = count
     
  seq_ids = sorted(seq_ids.keys())
  
  info('Found %d sequences at least %d kb in length' % (len(seq_ids), min_kb_size))
  
  name_dict = get_chromosome_names(seq_ids, email)
  
  with open(out_path, 'w') as file_obj:
    for seq_id in sorted(name_dict):
      line = '%s\t%s\n' % (seq_id, name_dict[seq_id])
      file_obj.write(line)
  
  info('Written %d chromsome names to "%s"' % (len(name_dict), out_path))


def main(argv=None):

  from argparse import ArgumentParser

  if argv is None:
    argv = sys.argv[1:]

  epilog = 'For further help email tjs23@cam.ac.uk or wb104@cam.ac.uk'

  arg_parse = ArgumentParser(prog=PROG_NAME, description=DESCRIPTION,
                             epilog=epilog, prefix_chars='-', add_help=True)
   
  arg_parse.add_argument(metavar='OUT_FILE', dest='o',
                         help='The TSV format output file that chromosome names will be written to')
 
  arg_parse.add_argument(metavar='EMAIL', dest='e',
                         help='An email address to tag remote NCBI Entrez lookups')

  arg_parse.add_argument(nargs='+', metavar='FASTA_FILE', dest='f',
                         help='Genome FASTA files, exactly as used by nuc_process, for a genome build' \
                              ' (accepts wildcards)')

  arg_parse.add_argument('-m', '--min-size', default=DEFAULT_MIN_SIZE, metavar='MIN_CHR_SIZE', dest='m',
                         type=int, help='Minimum sequence length to consider in kb. Default: %d' % DEFAULT_MIN_SIZE)
  
  args = vars(arg_parse.parse_args(argv))
  
  email = args['e']
  out_path = args['o']
  fasta_paths = args['f']
  min_size = args['m']
  
  write_chromosome_name_file(out_path, email, fasta_paths, min_size)
  

if __name__ == "__main__":
  main()
