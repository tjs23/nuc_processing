import sys
import os
import gzip

from collections import defaultdict
from sys import stdout

PROG_NAME = 'split_fastq_barcodes'
VERSION = '1.0.0'
DESCRIPTION = 'A Python script to split paired read FASTAQ files into separate files according to thier 3\' and 5\' barcode sequences'
IO_BUFFER = int(1e6)
DEFAULT_MIN_COUNT = 1000
DEFAULT_BC_LEN = 3
COUNT_LEVELS = (1000000, 100000, 10000, 1000, 100)
 
def check_regular_file(file_path):

  msg = ''
  
  if not os.path.exists(file_path):
    msg = 'File "%s" does not exist'
    return False, msg % file_path
  
  if not os.path.isfile(file_path):
    msg = 'Location "%s" is not a regular file'
    return False, msg % file_path
  
  if os.stat(file_path).st_size == 0:
    msg = 'File "%s" is of zero size '
    return False, msg % file_path
    
  if not os.access(file_path, os.R_OK):
    msg = 'File "%s" is not readable'
    return False, msg % file_path
  
  return True, msg
  

def report(msg):

  print(msg)
   

def warn(msg, prefix='WARNING'):

  report('%s: %s' % (prefix, msg))

 
def critical(msg, prefix='FAILURE'):

  report('%s: %s' % (prefix, msg))
  sys.exit(0)


def info(msg, prefix='INFO'):

  report('%s: %s' % (prefix, msg))


def open_file(file_path, mode=None, gzip_exts=('.gz','.gzip')):
  """
  GZIP agnostic file opening
  """
  
  if os.path.splitext(file_path)[1].lower() in gzip_exts:
    file_obj = gzip.open(file_path, mode or 'rt')
  else:
    file_obj = open(file_path, mode or 'rU', IO_BUFFER)
  
  return file_obj
  

def split_fastq_barcodes(fastq_paths, bc_file_path=None, bc_length=3, out_dir=None, diff_ends=False, min_seqs=DEFAULT_MIN_COUNT, buff_size=1000):
             
  if len(fastq_paths) != 2:
    msg = 'Exectly two FASTQ files required (%d specified)' % len(fastq_paths)
    critical(msg)
  
  if out_dir:
    if not os.path.exists(out_dir):
      msg = 'Specified output directory "%s" does not exist' % out_dir
      critical(msg)

  else:
    out_dir = os.path.dirname(fastq_paths[0])
  
  for file_path in fastq_paths:
    is_ok, msg = check_regular_file(file_path)
    
    if not is_ok:
      critical(msg)
  
  sample_names = {}
  
  if bc_file_path:
    bc_length = None
    
    is_ok, msg = check_regular_file(bc_file_path)
    
    if not is_ok:
      critical(msg)
  
    with open(bc_file_path) as file_obj:
      for line in file_obj:
        vals = line.split()
        
        if not vals:
          continue
        
        if vals[0][0] == '#':
          continue
        
        bc_key = '_'.join(vals[:-1])
        
        #if bc_key in sample_names:
        #  msg = 'Barcode "%s" repeated' % ' '.join(vals[:-1])
        #  critical(msg)   
               
        sample_names[bc_key] = vals[-1]
  
      # Look at last one
      bc_length = len(vals[0])
      
      # Check consistencty of length single/dual
      if len(vals) > 2:
        diff_ends = True

        for bc_key in sample_names:
          if '_' not in bc_key:
            msg = 'Barcodes not consistently dual or single'
            critical(msg)
          
          if len(bc_key) != 2*bc_length + 1:
            msg = 'Barcodes not consistent length'
            critical(msg)

        info('Found %d barcode pairs in file %s' % (len(sample_names), bc_file_path))
      
      else:
        diff_ends = False
        
        for bc_key in sample_names:
          if '_' in bc_key:
            msg = 'Barcodes not consistently dual or single'
            critical(msg)
          
          if len(bc_key) != bc_length:
            msg = 'Barcodes not consistent length'
            critical(msg)
    
        info('Found %d barcodes in file %s' % (len(sample_names), bc_file_path))
     
      info('Barcode length %d' % bc_length)
     
  buff_size = max(min_seqs, buff_size)
  
  fastq_path_1, fastq_path_2 = fastq_paths
  
  file_name_1 = os.path.basename(fastq_path_1)
  file_name_2 = os.path.basename(fastq_path_2)
  
  if file_name_1.endswith('.gz'):
    file_name_1 = file_name_1[:-3]

  if file_name_2.endswith('.gz'):
    file_name_2 = file_name_2[:-3]  
    
  file_root_1, file_ext_1 = os.path.splitext(file_name_1)
  path_root_1 = os.path.join(out_dir, file_root_1)

  file_root_2, file_ext_2 = os.path.splitext(file_name_2)
  path_root_2 = os.path.join(out_dir, file_root_2)

  out_data_1 = defaultdict(list)
  out_data_2 = defaultdict(list)
  bc_counts = defaultdict(int)
  start_idx = bc_length + 1
  n_buff = defaultdict(int)
  n_pairs = 0
  n_valid = 0
  n_mismatched = 0
  level_counts = {x:0 for x in COUNT_LEVELS}
  
  def get_bc_file_name(path_root, bc_key, sample_names, file_ext):
    return '%s_bc_%s%s' % (path_root, sample_names.get(bc_key, bc_key), file_ext)
  
  with open_file(fastq_path_1) as file_obj_1, open_file(fastq_path_2) as file_obj_2:
    
    readline_1 = file_obj_1.readline
    readline_2 = file_obj_2.readline

    nn1 = 0
    nn2 = 0
    
    line_1a = readline_1()
    line_1b = readline_1()
    line_1c = readline_1()
    line_1d = readline_1()    
    line_2a = readline_2()
    line_2b = readline_2()
    line_2c = readline_2()
    line_2d = readline_2()    
 
    while line_1a:
      if line_1a[0] != '@':
        msg = 'FASTQ file "%s" line does not start with @: "%s"' % (fastq_path_1, line_1a)
        critical(msg)

      if line_2a[0] != '@':
        msg = 'FASTQ file "%s" line does not start with @: "%s"' % (fastq_path_2, line_2a)
        critical(msg)
          
      if line_1b[bc_length] == 'N':
        nn1 += 1
 
      if line_2b[bc_length] == 'N':
        nn2 += 1
 
      barcode_1 = line_1b[:bc_length]
      barcode_2 = line_2b[:bc_length]
      n_pairs += 1
      
      if barcode_1 != barcode_2 and not diff_ends:
        n_mismatched += 1

      else:
      
        if diff_ends:
          bc_key = '%s_%s' % (barcode_1, barcode_2)
        
        else:
          bc_key = barcode_1
       
        bc_counts[bc_key] += 1
        n = bc_counts[bc_key]
        
        if not sample_names or bc_key in sample_names:
          n_valid += 1
          out_data_1[bc_key] += [line_1a+line_1b[start_idx:]+line_1c+line_1d[start_idx:]]
          out_data_2[bc_key] += [line_2a+line_2b[start_idx:]+line_2c+line_2d[start_idx:]]
          n_buff[bc_key] += 1
          
          if n_buff[bc_key] > buff_size:
            if n < 2 * buff_size: # First time
              mode = 'w'
            else:
              mode = 'a'
 
            # Number of barcodes may exceed max number of open files
 
            out_file_obj_1 = open(get_bc_file_name(path_root_1, bc_key, sample_names, file_ext_1), mode, IO_BUFFER)
            out_file_obj_2 = open(get_bc_file_name(path_root_2, bc_key, sample_names, file_ext_2), mode, IO_BUFFER)
 
            for line4 in out_data_1[bc_key]:
              out_file_obj_1.write(line4)
 
            for line4 in out_data_2[bc_key]:
              out_file_obj_2.write(line4)
 
            out_file_obj_1.close()
            out_file_obj_2.close()

            out_data_1[bc_key] = []
            out_data_2[bc_key] = []
            n_buff[bc_key] = 0
        
        if n in level_counts:
          level_counts[n] += 1
        
      line_1a = readline_1()
      line_1b = readline_1()
      line_1c = readline_1()
      line_1d = readline_1()
      line_2a = readline_2()
      line_2b = readline_2()
      line_2c = readline_2()
      line_2d = readline_2()
      
      if n_pairs % 100000 == 0:
        levels = ' '.join(['{:,}>{:,}'.format(level_counts[l], l) for l in COUNT_LEVELS])
        n_bc  = len(bc_counts)
        stdout.write("\r  Reads:{:,} Valid:{:,} BCs:{:,} - {}".format(n_pairs, n_valid, n_bc, levels))
        stdout.flush()
    
    for bc_key in bc_counts:
      if not sample_names and bc_counts[bc_key] < min_seqs: # Ignore small counts if not using a preset list of barcodes
        continue
      
      if not n_buff[bc_key]:
        continue
      
      if bc_counts[bc_key] < 2 * buff_size: # First time
        mode = 'w'
      else:
        mode = 'a'
      
      # Write anything that remains
      
      out_file_obj_1 = open(get_bc_file_name(path_root_1, bc_key, sample_names, file_ext_1), mode, IO_BUFFER)
      out_file_obj_2 = open(get_bc_file_name(path_root_2, bc_key, sample_names, file_ext_2), mode, IO_BUFFER)
   
      for line4 in out_data_1[bc_key]:
        out_file_obj_1.write(line4)
         
      for line4 in out_data_2[bc_key]:
        out_file_obj_2.write(line4)
      
      out_file_obj_1.close()
      out_file_obj_2.close()
    
      out_data_1[bc_key] = []
      out_data_2[bc_key] = []
      n_buff[bc_key] = 0

  stdout.write("\n") # move the cursor to the next line
    
  if diff_ends:
    msg = 'Found %d read pairs' % n_pairs
  
  else:
    msg = 'Found %d read pairs and %d mismatches (%.3f%%)' % (n_pairs, n_mismatched, n_mismatched/float(n_pairs))
    
  info(msg)  
  
  msg = 'Found unknown base (N) at start of %d and %d reads for respective inputs' % (nn1, nn2)
  
  info(msg)  
  
  if diff_ends:
    for bc_key in sorted(bc_counts):
      nb = bc_counts[bc_key]
      msg = 'Barcode: {} count {:,} ({:.2f}%)'.format(bc_key, nb, nb/float(100.0*n_pairs))
      info(msg)  

  else:
    for bc_key in sorted(bc_counts):
      nb = bc_counts[bc_key]
      bc1, bc2 = bc_key.split('_')
      msg = 'Barcodes: {} {} count {:,} ({:.2f}%)'.format(bc1, bc2, nb, nb/float(100.0*n_pairs))
      info(msg)  


if __name__ == '__main__':
 
  from argparse import ArgumentParser

  epilog = 'For further help email tjs23@cam.ac.uk or wb104@cam.ac.uk'
 
  arg_parse = ArgumentParser(prog=PROG_NAME, description=DESCRIPTION,
                            epilog=epilog, prefix_chars='-', add_help=True)  

  arg_parse.add_argument('i', nargs=2, metavar='TWO_FASTQ_FILES',
                         help='The two input paired-read FASTQ files to process. Accepts wildcards that match a pair of files')

  arg_parse.add_argument('-bc', metavar='BARCODE_SAMPLE_TSV', default=None,
                         help='Optional TSV file containing vaid barcodes in first column (or first two columns in R1, R2 order) and corresponding sample names for output in last column.')

  arg_parse.add_argument('-o', metavar='OUT_DIR', default=None,
                         help='Output directory for the results. Defaults to the directory for the first input FASTQ file.')

  arg_parse.add_argument('-m',  '--min-seqs', metavar='MIN_COUNT', dest='m', default=DEFAULT_MIN_COUNT, type=int,
                         help='Minimum number of sequence entries required before a barcoded set is written to file. Ignored if using "bc" option. Default: %d' % DEFAULT_MIN_COUNT)

  arg_parse.add_argument('-s',  '--barcode-size', metavar='BARCODE_SIZE', dest='s', default=DEFAULT_BC_LEN, type=int,
                         help='Barcode length in basepairs. Ignored if using "bc" option. Default: %d'% DEFAULT_BC_LEN)

  arg_parse.add_argument('-d', '--different-ends', action='store_true', dest='d',
                         help='Specifies that potentially different barcodes are used for the 5\' and 3\' ends of the same read pair. Ignored if using "bc" option.')    
  
  args = vars(arg_parse.parse_args())
  
  fastq_paths  = args['i']
  bc_length    = args['s']
  out_dir      = args['o']
  diff_ends    = args['d']
  min_seqs     = args['m']
  bc_file_path = args['bc']
  
  split_fastq_barcodes(fastq_paths, bc_file_path, bc_length, out_dir, diff_ends, min_seqs)
