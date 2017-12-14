"""
---- COPYRIGHT ----------------------------------------------------------------

Copyright (C) 20016-2017
Tim Stevens (MRC-LMB) and Wayne Boucher (University of Cambridge)


---- LICENSE ------------------------------------------------------------------

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation, either version 3 of the License, or (at your option) any
later version.

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

import sys
import os


def demultiplex(barcode_file, fastq_paths, buff_size=10000):
  """
  Script to separate single-cell data from the combined massively multiplexed Hi-C
  data described by:
  Massively multiplex single-cell Hi-C. Nature Methods 14, 263-3266 (2017)
  """
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
This script is run after doing the equivalent to:

 SeqPrep/SeqPrep -A AGATCGGAAGAGCGATCGG -B AGATCGGAAGAGCGTCGTG -f SRR3956928_1.fastq.gz -r SRR3956928_2.fastq.gz -1 SRR3956928_1.fastq.clipped -2 SRR3956928_2.fastq.clipped > seq_prep.txt 2>> adaptor_clipping_stats.txt

 python inline_splitter.py ../SRR3956928_1_clip.fastq.gz ../SRR3956928_2_clip.fastq.gz outer_barcodes.txt ../SRR3956928_1_split ../SRR3956928_2_split

 python analyze_scDHC_V2design.py inner_barcodes.txt ../SRR3956928_1_split ../SRR3956928_2_split ../SRR3956928_final_r1.fastq ../SRR3956928_final_r2.fastq > ../SRR3956928_demultiplex

as available at https://github.com/VRam142/combinatorialHiC
"""
