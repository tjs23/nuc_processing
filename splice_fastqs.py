from sys import stdout, stderr, exit
from time import time

PROG_NAME = 'splice_fastqs'
VERSION = '1.0.0'
DESCRIPTION = 'A Python script to join sequence and quality data (e.g. barcodes to genomic sequence) in two ordered FASTQ files'


def splice_fastqs(fastq_path_a, fastq_path_b, out_path=None, check_ids=True, keep_rep_ids=False, io_buffer=int(1e8)): # 2,147,483,647
  
  if out_path:
    out_file_obj = open(out_path, 'w', buffering=io_buffer)
  else:
    out_file_obj = stdout
  
  write = out_file_obj.write
  i = 0
  
  t0 = time()
    
  with open(fastq_path_a, 'r', buffering=io_buffer) as fq_a, open(fastq_path_b, buffering=io_buffer) as fq_b:
    
    out_lines = []
    
    read_a = fq_a.readline
    read_b = fq_b.readline
    
    line_a1 = read_a()
    line_a2 = read_a()
    line_a3 = read_a()
    line_a4 = read_a()
    line_b1 = read_b()
    line_b2 = read_b()
    line_b3 = read_b()
    line_b4 = read_b()
    
    while line_a1 and line_b1:
      
      if check_ids:
        m = line_a1.rfind(' ')
        n = line_b1.rfind(' ')
        if line_a1[:m] != line_b1[:n]:
          stderr.write('FAILURE: Line mismatch at entry %d when comparing %s and %s\n' % (i, fastq_path_a, fastq_path_b))
          stderr.write(line_a1)
          stderr.write(line_b1)
          exit(0)
      
      line_1 = line_a1
      line_2 = line_a2[:-1] + line_b2
      
      if keep_rep_ids:
        line_3 = line_a3
      
      else:
        line_3 = '+\n'
      
      line_4 = line_a4[:-1] + line_b4
      
      write(line_1)
      write(line_2)
      write(line_3)
      write(line_4)
      
      i += 1
      if i % 100000 == 0:
        if out_path:
          stdout.write("\r  Processed {:,} entries in {:.1f}s".format(i, time()-t0))
          stdout.flush()
      
      line_a1 = read_a()
      line_a2 = read_a()
      line_a3 = read_a()
      line_a4 = read_a()
      line_b1 = read_b()
      line_b2 = read_b()
      line_b3 = read_b()
      line_b4 = read_b()
      
  if out_lines:
    write_lines(out_lines)

  if out_path:
    stdout.write("\r  Processed {:,} entries in {:.1f}s\n".format(i, time()-t0))
    stdout.flush()
    out_file_obj.close()
  

if __name__ == '__main__':
 
  from argparse import ArgumentParser

  epilog = 'For further help email tjs23@cam.ac.uk or wb104@cam.ac.uk'
 
  arg_parse = ArgumentParser(prog=PROG_NAME, description=DESCRIPTION,
                             epilog=epilog, prefix_chars='-', add_help=True)  

  arg_parse.add_argument('fastqs', nargs=2, metavar='FASTQ_FILES',
                         help='Two input FASTQ of equal length and with corresponding entry order. The sequences from the first file will be joined to the beginning of sequences in the second file.')
  
  arg_parse.add_argument('-o', metavar='OUT_FILE', default=None,
                         help='Optional name for output FASTQ file. By default output will be printed to STDOUT')

  arg_parse.add_argument('-c', '--check-ids', action='store_true', dest='c',
                         help="Check that joined FASTQ entries have the same read IDs. When not set the IDs from the first FASTQ file will be used.")    

  arg_parse.add_argument('-k', '--keep-repeat-ids', action='store_true', dest='k',
                         help='Specified that any repeat sequence IDs (on FASTQ lines starting with "+") should be kept. Keeping repeat IDs will result in a larger file')
  
  args = vars(arg_parse.parse_args())
  
  fastq_path_a, fastq_path_b  = args['fastqs']
  
  out_path     = args['o']
  check_ids    = args['c']
  keep_rep_ids = args['k']
  
  splice_fastqs(fastq_path_a, fastq_path_b, out_path, check_ids, keep_rep_ids)
