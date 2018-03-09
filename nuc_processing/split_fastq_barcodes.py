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
import gzip

from collections import defaultdict
from sys import stdout

PROG_NAME = 'split_fastq_barcodes'
VERSION = '1.0.0'
DESCRIPTION = 'A Python script to split FASTQ files, both paired or single-end, into separate files according to thier 5\' and 3\' barcode sequences'
IO_BUFFER = int(4e6)
ANALYSIS_TAG = 'bc_report'
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


def read_barcode_file(bc_file_path):

  sample_names = {}

  sample = 0

  with open(bc_file_path) as file_obj:
    for line in file_obj:
      vals = line.split()

      if not vals:
        continue

      if vals[0][0] == '#':
        continue

      bc_key = vals[0]

      if bc_key in sample_names:
        msg = 'Barcode "%s" repeated' % bc_key
        critical(msg)

      if len(vals) > 1:
        sample_names[bc_key] = vals[1]

      else:
        sample += 1
        sample_names[bc_key] = 'sample_%d' % sample

    # Look at last one
    bc_length = len(bc_key)

    # Check consistencty of length single/dual
    if '-' in bc_key:
      for bc_key2 in sample_names:
        if '-' not in bc_key2:
          msg = 'Barcodes not consistently dual or single'
          critical(msg)

        if len(bc_key) != bc_length:
          msg = 'Barcodes not consistent length'
          critical(msg)

      info('Found %d barcode pairs in file %s' % (len(sample_names), bc_file_path))
      info('Barcode length %d' % ((bc_length-1)/2))

    else:
      for bc_key2 in sample_names:
        if '-' in bc_key2:
          msg = 'Barcodes not consistently dual or single'
          critical(msg)

        if len(bc_key) != bc_length:
          msg = 'Barcodes not consistent length'
          critical(msg)

      info('Found %d barcodes in file %s' % (len(sample_names), bc_file_path))
      info('Barcode length %d' % bc_length)

  return sample_names, len(bc_key.split('-')[0]), '-' in bc_key


def _write_analysis_file(bc_counts, sample_names, analysis_file_path):

  data_list = []
  n = 0

  for bc in bc_counts:
    count = bc_counts[bc]
    sample = sample_names.get(bc)

    n += count
    data_list.append([count, bc, sample])

  data_list.sort(reverse=True)
  n = float(n)
  s = 0

  with open(analysis_file_path, 'w') as file_obj:
    line = '#%s\t%s\t%s\t%s\n' % ('barcode', 'sample_name', 'count', '%abundance')
    file_obj.write(line)

    for count, bc, sample in data_list:
      if sample is None:
        s += 1
        sample = 'sample_%d' % s

      line = '{}\t{}\t{:,}\t{:.2f}\n'.format(bc, sample, count, 100.0 * count/n)
      file_obj.write(line)

  info('Barcode analysis written to %s' % analysis_file_path)


def analyse_fastq_barcodes(fastq_paths, sample_names, analysis_file_path, bc_length, diff_ends, in_illumina_head):

  n_reads = 0
  bc_counts = defaultdict(int)
  is_paired = len(fastq_paths) == 2
  level_counts = {x:0 for x in COUNT_LEVELS}

  if is_paired:
    fastq_path_1, fastq_path_2 = fastq_paths

    with open_file(fastq_path_1) as file_obj_1, open_file(fastq_path_2) as file_obj_2:
      readline_1 = file_obj_1.readline
      readline_2 = file_obj_2.readline

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

        if in_illumina_head:
          barcode_1, barcode_2 = line_1a.strip().split(':')[-1].split('+')

        else:
          barcode_1 = line_1b[:bc_length]
          barcode_2 = line_2b[:bc_length]

        n_reads += 1

        if barcode_1 != barcode_2 and not diff_ends:
          continue

        else:
          if diff_ends:
            bc_key = '%s-%s' % (barcode_1, barcode_2)
          else:
            bc_key = barcode_1

          bc_counts[bc_key] += 1
          n = bc_counts[bc_key]

          if n in level_counts:
            level_counts[n] += 1

        line_1a = readline_1()
        line_1b = readline_1()
        line_1c = readline_1()
        line_1d = readline_1()
        line_2a = readline_2()
        line_2b = readline_2()
        line_2c = readline_2() # noqa: F841
        line_2d = readline_2() # noqa: F841

        if n_reads % 100000 == 0:
          levels = ' '.join(['{:,}>{:,}'.format(level_counts[l], l) for l in COUNT_LEVELS])
          n_bc = len(bc_counts)
          stdout.write("\r  Reads:{:,} BC seqs:{:,} - {}".format(n_reads, n_bc, levels))
          stdout.flush()

  else: # Single-end reads: one file
    fastq_path_1 = fastq_paths[0]

    with open_file(fastq_path_1) as file_obj_1:
      readline_1 = file_obj_1.readline

      line_1a = readline_1()
      line_1b = readline_1()
      line_1c = readline_1()
      line_1d = readline_1()

      while line_1a:
        if line_1a[0] != '@':
          msg = 'FASTQ file "%s" line does not start with @: "%s"' % (fastq_path_1, line_1a)
          critical(msg)

        if in_illumina_head:
          barcode_1, barcode_2 = line_1a.strip().split(':')[-1].split('+')
          bc_key = '%s-%s' % (barcode_1, barcode_2)

        else:
          bc_key = line_1b[:bc_length]

        n_reads += 1
        bc_counts[bc_key] += 1
        n = bc_counts[bc_key]

        if n in level_counts:
          level_counts[n] += 1

        line_1a = readline_1()
        line_1b = readline_1()
        line_1c = readline_1() # noqa: F841
        line_1d = readline_1() # noqa: F841

        if n_reads % 100000 == 0:
          levels = ' '.join(['{:,}>{:,}'.format(level_counts[l], l) for l in COUNT_LEVELS])
          n_bc = len(bc_counts)
          stdout.write("\r  Reads:{:,} BC seqs:{:,} - {}".format(n_reads, n_bc, levels))
          stdout.flush()

  stdout.write("\n") # move the cursor to the next line

  _write_analysis_file(bc_counts, sample_names, analysis_file_path)


def split_fastq_barcodes(fastq_paths, bc_file_path=None, analysis_file_path=None, out_dir=None, max_mismatches=0,
                         bc_length=None, diff_ends=False, in_illumina_head=False, no_samp_name=False, write_unmatched=True,
                         buff_size=10000):

  if len(fastq_paths) not in (1,2):
    critical('One FASTQ file or two (paired) FASTQ files are required (%d specified)' % len(fastq_paths))

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

    if not analysis_file_path:
      file_name = os.path.basename(file_path)
      analysis_file_path = os.path.join(out_dir, '%s_%s.tsv' % (os.path.splitext(file_name)[0], ANALYSIS_TAG))

  is_paired = len(fastq_paths) == 2

  if bc_file_path:
    bc_length = None
    is_ok, msg = check_regular_file(bc_file_path)

    if not is_ok:
      critical(msg)

    sample_names, bc_length, diff_ends = read_barcode_file(bc_file_path)

  else:
    if not (bc_length or in_illumina_head):
      msg = 'Barcode length (-s) must be supplied when -b (barcode file) or -ih options are not specified.'
      critical(msg)

    info('No -b barcode file specified.')
    info('No FASTQ files will be written.')
    info('Only analysis will be performed.')
    analyse_fastq_barcodes(fastq_paths, {}, analysis_file_path, bc_length, diff_ends, in_illumina_head)

    return

  def get_bc_file_name(path_root, bc_key, sample_names, nsn, file_ext):
    if nsn:
      sn = bc_key
    else:
      sn = sample_names.get(bc_key, bc_key)

    return '%s_%s%s' % (path_root, sn, file_ext)

  bc_counts = defaultdict(int)
  n_buff = defaultdict(int)
  n_reads = 0
  n_valid = 0
  n_mismatched = 0
  n_lost = 0
  n_imperfect = 0
  level_counts = {x:0 for x in COUNT_LEVELS}

  if diff_ends:
    valid_1, valid_2 = zip(*[x.split('-') for x in sample_names])

  else:
    valid_1 = valid_2 = list(sample_names.keys())

  valid_1 = set(valid_1)
  valid_2 = set(valid_2)
  opened_files = set()

  if in_illumina_head:
    start_idx = 0
  else:
    start_idx = bc_length

  if is_paired:
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

    if write_unmatched:
      lost_path_1 = '%s_lost_reads%s' % (path_root_1, file_ext_1)
      lost_path_2 = '%s_lost_reads%s' % (path_root_2, file_ext_2)

      lost_file_obj_1 = open_file(lost_path_1, 'w')
      lost_file_obj_2 = open_file(lost_path_2, 'w')

    out_data_1 = defaultdict(list)
    out_data_2 = defaultdict(list)

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

        if line_1b[start_idx] == 'N':
          nn1 += 1

        if line_2b[start_idx] == 'N':
          nn2 += 1

        if in_illumina_head:
          barcode_1, barcode_2 = line_1a.strip().split(':')[-1].split('+')

        else:
          barcode_1 = line_1b[:bc_length]
          barcode_2 = line_2b[:bc_length]

        n_reads += 1

        if barcode_1 != barcode_2 and not diff_ends:
          n_mismatched += 1

        else:

          if diff_ends:
            bc_key = '%s-%s' % (barcode_1, barcode_2)
          else:
            bc_key = barcode_1

          bc_counts[bc_key] += 1
          n = bc_counts[bc_key]

          if bc_key in sample_names:
            valid = True

          elif max_mismatches:
            valid = False

            if barcode_1 in valid_1:
              poss_1 = [barcode_1]

            else:
              poss_1 = []
              for bc in valid_1:
                m = 0

                for i in range(bc_length):
                  if bc[i] != barcode_1[i]:
                    m += 1

                if m <= max_mismatches:
                  poss_1.append(bc)

            if len(poss_1) == 1:
              if diff_ends:
                if barcode_2 in valid_2:
                  poss_2 = [barcode_2]

                else:
                  poss_2 = []
                  for bc in valid_2:
                    m = 0

                    for i in range(bc_length):
                      if bc[i] != barcode_2[i]:
                        m += 1

                    if m <= max_mismatches:
                      poss_2.append(bc)

                if len(poss_2) == 1:
                  bc_key = '%s-%s' % (poss_1[0], poss_2[0])
                  n_imperfect += 1
                  valid = True

              else:
                bc_key = poss_1[0]
                n_imperfect += 1
                valid = True

          else:
            valid = False

          if valid:
            out_data_1[bc_key] += [line_1a+line_1b[start_idx:]+line_1c+line_1d[start_idx:]]
            out_data_2[bc_key] += [line_2a+line_2b[start_idx:]+line_2c+line_2d[start_idx:]]
            n_buff[bc_key] += 1
            n_valid += 1

            if n_buff[bc_key] > buff_size:

              file_path_1 = get_bc_file_name(path_root_1, bc_key, sample_names, no_samp_name, file_ext_1)
              file_path_2 = get_bc_file_name(path_root_2, bc_key, sample_names, no_samp_name, file_ext_2)

              if file_path_1 in opened_files:
                mode = 'a'
              else:
                opened_files.add(file_path_1)
                opened_files.add(file_path_2)
                mode = 'w'

              # Number of (valid) barcodes may exceed max number of open files
              out_file_obj_1 = open(file_path_1, mode, IO_BUFFER)
              out_file_obj_2 = open(file_path_2, mode, IO_BUFFER)

              for line4 in out_data_1[bc_key]:
                out_file_obj_1.write(line4)

              for line4 in out_data_2[bc_key]:
                out_file_obj_2.write(line4)

              out_file_obj_1.close()
              out_file_obj_2.close()

              out_data_1[bc_key] = []
              out_data_2[bc_key] = []
              n_buff[bc_key] = 0

          else: # 'Lost'
            if write_unmatched:
              lost_file_obj_1.write(line_1a+line_1b+line_1c+line_1d)
              lost_file_obj_2.write(line_2a+line_2b+line_2c+line_2d)

            n_lost += 1

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

        if n_reads % 100000 == 0:
          levels = ' '.join(['{:,}>{:,}'.format(level_counts[l], l) for l in COUNT_LEVELS])
          n_bc = len(bc_counts)

          if max_mismatches:
            stdout.write("\r  Reads:{:,} Perfect:{:,} Imperfect:{:,} BCs:{:,} - {}".format(n_reads, n_valid, n_imperfect, n_bc, levels))

          else:
            stdout.write("\r  Reads:{:,} Valid:{:,} BCs:{:,} - {}".format(n_reads, n_valid, n_bc, levels))

          stdout.flush()

      for bc_key in bc_counts: # Handle remaining cache
        if n_buff[bc_key]:
          file_path_1 = get_bc_file_name(path_root_1, bc_key, sample_names, no_samp_name, file_ext_1)
          file_path_2 = get_bc_file_name(path_root_2, bc_key, sample_names, no_samp_name, file_ext_2)

          if file_path_1 in opened_files:
            mode = 'a'
          else:
            mode = 'w'

          # Number of (valid) barcodes may exceed max number of open files
          out_file_obj_1 = open(file_path_1, mode, IO_BUFFER)
          out_file_obj_2 = open(file_path_2, mode, IO_BUFFER)

          for line4 in out_data_1[bc_key]:
            out_file_obj_1.write(line4)

          for line4 in out_data_2[bc_key]:
            out_file_obj_2.write(line4)

          out_file_obj_1.close()
          out_file_obj_2.close()
          out_data_1[bc_key] = []
          out_data_2[bc_key] = []
          n_buff[bc_key] = 0

    if write_unmatched:
      lost_file_obj_1.close()
      lost_file_obj_2.close()

  else: # Single-end reads: one file
    fastq_path_1 = fastq_paths[0]
    file_name_1 = os.path.basename(fastq_path_1)

    if file_name_1.endswith('.gz'):
      file_name_1 = file_name_1[:-3]

    file_root_1, file_ext_1 = os.path.splitext(file_name_1)
    path_root_1 = os.path.join(out_dir, file_root_1)
    out_data_1 = defaultdict(list)

    if write_unmatched:
      lost_path_1 = '%s_lost_reads%s' % (path_root_1, file_ext_1)
      lost_file_obj_1 = open_file(lost_path_1, 'w')

    with open_file(fastq_path_1) as file_obj_1:
      readline_1 = file_obj_1.readline
      nn1 = 0

      line_1a = readline_1()
      line_1b = readline_1()
      line_1c = readline_1()
      line_1d = readline_1()

      while line_1a:
        if line_1a[0] != '@':
          msg = 'FASTQ file "%s" line does not start with @: "%s"' % (fastq_path_1, line_1a)
          critical(msg)

        if line_1b[start_idx] == 'N':
          nn1 += 1

        if in_illumina_head:
          barcode_1, barcode_2 = line_1a.strip().split(':')[-1].split('+')
          bc_key = '%s-%s' % (barcode_1, barcode_2)

        else:
          bc_key = barcode_1 = line_1b[:bc_length]

        n_reads += 1
        bc_counts[bc_key] += 1
        n = bc_counts[bc_key]

        if bc_key in sample_names:
          valid = True

        elif max_mismatches:
          valid = False

          if barcode_1 in valid_1:
            poss_1 = [barcode_1]

          else:
            poss_1 = []
            for bc in valid_1:
              m = 0

              for i in range(bc_length):
                if bc[i] != barcode_1[i]:
                  m += 1

              if m <= max_mismatches:
                poss_1.append(bc)

          if len(poss_1) == 1:
            if in_illumina_head:
              if barcode_2 in valid_2:
                poss_2 = [barcode_2]

              else:
                poss_2 = []
                for bc in valid_2:
                  m = 0

                  for i in range(bc_length):
                    if bc[i] != barcode_2[i]:
                      m += 1

                  if m <= max_mismatches:
                    poss_2.append(bc)

              if len(poss_2) == 1:
                bc_key = '%s-%s' % (poss_1[0], poss_2[0])
                n_imperfect += 1
                valid = True

            else:
              bc_key = poss_1[0]
              n_imperfect += 1
              valid = True

        else:
          valid = False

        if valid:
          out_data_1[bc_key] += [line_1a+line_1b[start_idx:]+line_1c+line_1d[start_idx:]]
          n_buff[bc_key] += 1

          if n_buff[bc_key] > buff_size:
            file_path_1 = get_bc_file_name(path_root_1, bc_key, sample_names, no_samp_name, file_ext_1)

            if file_path_1 in opened_files:
              mode = 'a'
            else:
              opened_files.add(file_path_1)
              mode = 'w'

            # Number of (valid) barcodes may exceed max number of open files
            out_file_obj_1 = open(file_path_1, mode, IO_BUFFER)
            for line4 in out_data_1[bc_key]:
              out_file_obj_1.write(line4)

            out_file_obj_1.close()
            out_data_1[bc_key] = []
            n_buff[bc_key] = 0

        else: # Lost
          if write_unmatched:
            lost_file_obj_1.write(line_1a+line_1b+line_1c+line_1d)
          n_lost += 1

        if n in level_counts:
          level_counts[n] += 1

        line_1a = readline_1()
        line_1b = readline_1()
        line_1c = readline_1()
        line_1d = readline_1()

        if n_reads % 100000 == 0:
          levels = ' '.join(['{:,}>{:,}'.format(level_counts[l], l) for l in COUNT_LEVELS])
          n_bc = len(bc_counts)

          if max_mismatches:
            stdout.write("\r  Reads:{:,} Perfect:{:,} Imperfect:{:,} BCs:{:,} - {}".format(n_reads, n_valid, n_imperfect, n_bc, levels))
          else:
            stdout.write("\r  Reads:{:,} Valid:{:,} BCs:{:,} - {}".format(n_reads, n_valid, n_bc, levels))

          stdout.flush()

      for bc_key in bc_counts:
        if n_buff[bc_key]:
          file_path_1 = get_bc_file_name(path_root_1, bc_key, sample_names, no_samp_name, file_ext_1)

          if file_path_1 in opened_files:
            mode = 'a'
          else:
            mode = 'w'

          out_file_obj_1 = open(file_path_1, mode, IO_BUFFER)
          for line4 in out_data_1[bc_key]:
            out_file_obj_1.write(line4)

          out_file_obj_1.close()
          out_data_1[bc_key] = []
          n_buff[bc_key] = 0

    if write_unmatched:
      lost_file_obj_1.close()

  stdout.write("\n") # move the cursor to the next line

  if diff_ends:
    msg1 = 'Found %d read pairs' % n_reads
  else:
    msg1 = 'Found %d read pairs and %d end mismatches (%.3f%%)' % (n_reads, n_mismatched, n_mismatched/float(n_reads))

  if is_paired:
    msg2 = 'Found unknown base (N) at start of %d and %d reads for respective inputs' % (nn1, nn2)

  else:
    msg2 = 'Found unknown base (N) at start of %d reads' % (nn1, )

  levels = ' '.join(['{:,}>{:,}'.format(level_counts[l], l) for l in COUNT_LEVELS])
  msg3 = 'Barcode sequences: Expected:{:,} Found:{:,} - {}'.format(len(sample_names), len(bc_counts), levels)
  msg4 = 'Barcode counts: Perfect:{:,} Imperfect:{:,} Unallocated:{:,}'.format(n_valid, n_imperfect, n_lost)

  info(msg1)
  info(msg2)
  info(msg3)
  info(msg4)

  if diff_ends:
    for bc_key in sorted(sample_names):
      nb = bc_counts[bc_key]
      bc1, bc2 = bc_key.split('_')
      msg = '{} barcodes: {} {} count {:,} ({:.2f}%)'.format(sample_names[bc_key], bc1, bc2, nb, nb/float(100.0*n_reads))
      info(msg)

  else:
    for bc_key in sorted(sample_names):
      nb = bc_counts[bc_key]
      msg = '{} barcode: {} count {:,} ({:.2f}%)'.format(sample_names[bc_key], bc_key, nb, nb/float(100.0*n_reads))
      info(msg)

  _write_analysis_file(bc_counts, sample_names,  analysis_file_path)


if __name__ == '__main__':

  from argparse import ArgumentParser

  # -k Keep barcode seqs, don't trim

  # allow separate barcode files

  # # To consider
  # Seek complements -rc1 -rc2
  # GZIP output?
  # Use Illumina contents file?
  # Quality based truncation? Trim short?

  epilog = 'For further help email tjs23@cam.ac.uk or wb104@cam.ac.uk'

  arg_parse = ArgumentParser(prog=PROG_NAME, description=DESCRIPTION,
                             epilog=epilog, prefix_chars='-', add_help=True)

  arg_parse.add_argument('i', nargs='+', metavar='FASTQ_FILE(S)',
                         help='One or two input FASTQ files to process. If two files are input they are assumed to be paired reads with matching rows. Accepts wildcards that match two files')

  arg_parse.add_argument('-b', metavar='BARCODE_FILE', default=None,
                         help='Optional TSV file specifying expected/valid barcodes. Contains single barcodes or pairs joined with "-", one per row with optional sample names. An existing analysis output file (see -a option) may be used. If not specified only barcode analysis will be performed; split FASTQ files will not be created.')

  arg_parse.add_argument('-a', metavar='ANALYSIS_FILE', default=None,
                         help='Optional TSV file name for output of barcode analysis. Defaults to using a name derived from the input FASTQ files tagged with "%s"' % ANALYSIS_TAG)

  arg_parse.add_argument('-o', metavar='OUT_DIR', default=None,
                         help='Output directory for the results. Defaults to the directory for the first input FASTQ file.')

  arg_parse.add_argument('-m', '--max-mismatches', metavar='NUM_MISMATCHES', dest='m', default=0, type=int,
                         help='Maximum number of basepair mitchmatches tolerated in a barcode sequence, unless the mismatch is ambiguous (does not distinguish barcodes). Default: 0')

  arg_parse.add_argument('-s', '--barcode-size', metavar='BARCODE_SIZE', dest='s', default=None, type=int,
                         help='Barcode length in basepairs for analysis. Only required when not using -b or -ih options.')

  arg_parse.add_argument('-u', '--write-unmached', action='store_true', dest='u',
                         help='Write out FASTQ data which cannot be matched to a barcode (i.e. \'lost\' reads).')

  arg_parse.add_argument('-d', '--different-ends', action='store_true', dest='d',
                         help='Specifies that potentially different barcodes are used for the 5\' and 3\' ends of the same read pair. Ignored if using -b option.')

  arg_parse.add_argument('-nb', '--barcode-file-names', action='store_true', dest='nb',
                         help='Suppress the use of sample names in output FASTQ files and name files with barcode sequences instead.')

  arg_parse.add_argument('-ih', '--illimina-head-barcodes', action='store_true', dest='ih',
                         help='Take barcode sequences from Illumina format FASTQ header lines, rather than the sequence line. Assumes ":{barcode1}+{barcode2}" at end of the header.')

  args = vars(arg_parse.parse_args())

  fastq_paths = args['i']
  bc_length = args['s'] or None # Not zero
  out_dir = args['o']
  diff_ends = args['d']
  bc_file_path = args['b']
  an_file_path = args['a']
  in_ill_head = args['ih']
  no_samp_name = args['nb']
  max_mismatch = args['m']
  write_lost = args['u']

  split_fastq_barcodes(fastq_paths, bc_file_path, an_file_path, out_dir, max_mismatch, bc_length, diff_ends, in_ill_head, no_samp_name, write_lost)
