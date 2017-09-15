import sys
import os

def demultiplex(barcode_file, fastq_paths, buff_size=10000):
  
  print('Reading barcodes')
  
  barcode_dict = {}
  with open(barcode_file) as file_obj:
    for line in file_obj:
      tag, inner_a, inner_b, outer = line.split()
      
      if inner_a == inner_b:
        barcode_dict[tag] = (inner_a, outer)
    
  print('Found %d cell tags' % len(barcode_dict))
       
  for fastq_path in fastq_paths:
    print('Processing FASTQ file %s' % fastq_path)
    n = 0
    file_root, file_ext = os.path.splitext(fastq_path)
    out_data = {cell_key:[] for cell_key in barcode_dict.values()}
  
    with open(fastq_path) as file_obj:
      readline = file_obj.readline
      line1 = readline()
      
      while line1:
        line2 = readline()
        line3 = readline()
        line4 = readline()             
                
        tag = line1[1:].strip()
        cell_key = barcode_dict.get(tag)
        
        if cell_key:
          lines = (line1, line2, line3, line4)
          out_data[cell_key].append(lines)
 
          if len(out_data[cell_key]) == buff_size:
            inner, outer = cell_key
            out_path = '%s_%s_%s%s' % (file_root, inner, outer, file_ext)
            print('Dump %s' % out_path)
 
            with open(out_path, 'w') as out_file_obj:
 
              for lines in out_data[cell_key]:
                 out_file_obj.writelines(lines)

        n += 1
 
        if n % buff_size == 0:
          print('  .. %d' % n)
      
        line1 = readline()
        
    for cell_key in out_data:
      inner, outer = cell_key
      out_path = '%s_%s_%s%s' % (file_root, inner, outer, file_ext)
      print('Dump %s' % out_path)
          
      with open(out_path, 'w') as out_file_obj:
        for lines in out_data[cell_key]:
          if lines:
            out_file_obj.writelines(lines)
      
  
if __name__ == '__main__':

  paths = sys.argv[1:]

  barcode_file = paths[0]
  fastq_paths = paths[1:]
  
  demultiplex(barcode_file, fastq_paths)

"""
This scriptis run after doing the equivalent to:

 SeqPrep/SeqPrep -A AGATCGGAAGAGCGATCGG -B AGATCGGAAGAGCGTCGTG -f SRR3956928_1.fastq.gz -r SRR3956928_2.fastq.gz -1 SRR3956928_1.fastq.clipped -2 SRR3956928_2.fastq.clipped > seq_prep.txt 2>> adaptor_clipping_stats.txt

 python inline_splitter.py ../SRR3956928_1_clip.fastq.gz ../SRR3956928_2_clip.fastq.gz outer_barcodes.txt ../SRR3956928_1_split ../SRR3956928_2_split

 python analyze_scDHC_V2design.py inner_barcodes.txt ../SRR3956928_1_split ../SRR3956928_2_split ../SRR3956928_final_r1.fastq ../SRR3956928_final_r2.fastq > ../SRR3956928_demultiplex


as available at https://github.com/VRam142/combinatorialHiC

"""
