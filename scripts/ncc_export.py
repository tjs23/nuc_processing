
"""
---- COPYRIGHT ----------------------------------------------------------------

Copyright (C) 20016-2025
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
import sys, math, datetime, os, time, shutil
import numpy as np
import h5py

from collections import defaultdict
from nuc_tools.core import nuc_io as io
from nuc_tools.tools.ncc_bin import bin_ncc
from nuc_tools.core import nuc_util as util
from random import randint

PROG_NAME = 'ncc_export'
VERSION = '1.0.0'
DESCRIPTION = 'A Python script to export NCC format (RE fragment mapped ambiguous read level data) to common Hi-C data formats like .pairs, .cool etc.'
FORMATS = ('pairs','cool','mcool')
AVAIL_FORMATS = ', '.join(FORMATS)
DEFAULT_BIN_SIZE = 5.0
DEFAULT_BIN_SIZE_MCOOL = [1.0, 5.0, 10.0, 50.0, 100.0]
BUFFER_SIZE = 50000
REPORT_INTERVAL = 100_000
TEMP_FILE_FORMAT = '{}::{}::{}::{}.nccends'
TEMP_FILE_FORMAT_GZ = '{}::{}::{}::{}.nccends.gz'
CHUNK_SIZE = 50_000_000
 
def _preprocess_split_ncc(ncc_path, chunk_size=CHUNK_SIZE, line_buffer_size=BUFFER_SIZE):
    
    join = os.path.join
    
    print(f'INFO: Pre-read and split {ncc_path}')
    temp_dir = f'_ncc_export_temp_dir_{randint(1,10000000)}'
    os.mkdir(temp_dir)
    
    line_buffer = {}
    
    chromo_sizes = defaultdict(int)
    with io.open_file(ncc_path) as in_file_obj:
        for i, line in enumerate(in_file_obj):
            
            if i % REPORT_INTERVAL == 0:
                print(f'INFO: .. {i:,}', end='\r')
            
            chr_a, start_a, end_a, f_start_a, f_end_a, strand_a, \
              chr_b, start_b, end_b, f_start_b, f_end_b, strand_b, \
              ambig_group, pair_id, swap_pair = line.split()

            pos_a = int(f_start_a) if strand_a == '+' else int(f_end_a)
            pos_b = int(f_start_b) if strand_b == '+' else int(f_end_b)

            chromo_sizes[chr_a] = max(chromo_sizes[chr_a], pos_a)
            chromo_sizes[chr_b] = max(chromo_sizes[chr_b], pos_b)
 
            if chr_a > chr_b:
                chr_a, chr_b = chr_b, chr_a
                pos_a, pos_b = pos_b, pos_a
                strand_a, strand_b = strand_b, strand_a
            
            chunk_a = int(pos_a // chunk_size)
            chunk_b = int(pos_b // chunk_size)
            file_name = TEMP_FILE_FORMAT.format(chr_a, chr_b, chunk_a, chunk_b)
            
            if file_name not in line_buffer:
                line_buffer[file_name] = []
            
            line = f'{pair_id}\t{pos_a}\t{pos_b}\t{strand_a}\t{strand_b}\n'
            line_buffer[file_name].append(line)
            
            if len(line_buffer) >= line_buffer_size:
                file_path = join(temp_dir, file_name)
                with open(file_path, 'a') as out_file_obj:
                    out_file_obj.writelines(line_buffer[file_name])
               
                line_buffer[file_name] = []
                    
    for file_name in line_buffer:
        print(f'INFO: .. finalising {file_name}     ', end='\r')
        
        if line_buffer[file_name]:
            file_path = join(temp_dir, file_name)
            
            with open(file_path, 'a') as out_file_obj:
                out_file_obj.writelines(line_buffer[file_name])
            
            line_buffer[file_name] = []
            
        #io.compress_file(file_name)      
         
    print(f'INFO: .. {i:,}')
    
    return chromo_sizes, temp_dir
   

def _add_ncc_to_cool(ncc_path, hdf, bin_size, chromos, chromo_sizes, temp_dir):
    
    join = os.path.join
    
    hdf.attrs['format'] = 'HDF5::Cooler'
    hdf.attrs['format-version'] = 3
    hdf.attrs['bin-type'] = 'fixed'
    hdf.attrs['bin-size'] = bin_size
    hdf.attrs['storage-mode'] = 'symmetric-upper'
    hdf.attrs['generated-by'] = 'nuc_processing.ncc_export-{VERSION}'
    hdf.attrs['creation-date'] = datetime.datetime.now().isoformat()
    
    n_chromos = len(chromos)
    
    hdf_chroms = hdf.create_group('chroms')
    hdf_bins = hdf.create_group('bins')
    hdf_pixels = hdf.create_group('pixels')
    hdf_indexes = hdf.create_group('indexes')
    
    hdf_chroms_name = hdf_chroms.create_dataset('name', (n_chromos,), dtype='|S64', compression="gzip", data=np.array(chromos, dtype='S64'))
    hdf_chroms_length = hdf_chroms.create_dataset('length', (n_chromos,), dtype='i4', compression="gzip", data=[chromo_sizes[x] for x in chromos])
    
    chromo_offsets = {}
    chromo_idx = []
    chromo_starts = []
    chromo_ends = []
    offset = 0
    chromo_idx_bins = [offset]
    chromo_nbins = {}
    chromo_chunks = {}
    
    for k, chromo in enumerate(chromos):
        size = chromo_sizes[chromo]
        n_chr_bins = int(math.ceil(size/bin_size))
        chromo_nbins[chromo] = n_chr_bins 
        size_lim = bin_size * n_chr_bins
        starts = np.arange(0, size_lim, bin_size)
        ends = starts + (bin_size-1)
        chromo_chunks[chromo] = int(size/CHUNK_SIZE)
        
        chromo_idx.append(np.full(n_chr_bins, k))
        chromo_starts.append(starts)
        chromo_ends.append(ends)
        chromo_offsets[chromo] = offset
        offset += len(ends)
        chromo_idx_bins.append(offset)
        print(f'INFO: Chromosome {chromo:5s} size {size:12,d} regions {n_chr_bins:,}')
    
    chromo_idx = np.concatenate(chromo_idx, axis=0)
    chromo_starts = np.concatenate(chromo_starts, axis=0)
    chromo_ends = np.concatenate(chromo_ends, axis=0)
    n_bins = len(chromo_ends)
    
    hdf_bins.create_dataset('chrom', (n_bins,), dtype='i4', compression="gzip", data=chromo_idx)
    hdf_bins.create_dataset('start', (n_bins,), dtype='i4', compression="gzip", data=chromo_starts)
    hdf_bins.create_dataset('end',   (n_bins,), dtype='i4', compression="gzip", data=chromo_ends)

    #hdf_bins_weight = hdf_bins.create_dataset('weight', (n_bins,), dtype='i8', compression="gzip")
    
    n_pixels = n_bins * n_bins
    chromo_idx_bins = np.array(chromo_idx_bins, dtype=np.int64)
    bin_idx_pixels = np.arange(0, n_pixels+1, n_bins, dtype=np.int64)
    
    hdf_pixels_bin1_id = hdf_pixels.create_dataset('bin1_id', (n_pixels,), dtype='i8', compression="gzip", maxshape=(n_pixels,))
    hdf_pixels_bin2_id = hdf_pixels.create_dataset('bin2_id', (n_pixels,), dtype='i8', compression="gzip", maxshape=(n_pixels,))
    hdf_pixels_count   = hdf_pixels.create_dataset('count'  , (n_pixels,), dtype='i4', compression="gzip", maxshape=(n_pixels,))
    
    p1 = 0
    
    for k, chr_1 in enumerate(chromos):
        for chr_2 in chromos[k:]:
            if chr_1 > chr_2:
                 chr_a, chr_b = chr_2, chr_1
            else:
                 chr_a, chr_b = chr_1, chr_2
            
            t0 = time.time()
            print(f'INFO: Parsing {chr_1}:{chr_2}')
            t = 0
            
            for chunk_a in range(0, chromo_chunks[chr_a]+1):
                start_a = chunk_a * CHUNK_SIZE
                end_a = min(start_a + CHUNK_SIZE, chromo_sizes[chr_a])
                n_a = int(math.ceil((end_a-start_a)/bin_size)) # Num bins
                o_a = chromo_offsets[chr_a] + int(start_a/bin_size) # Bin offset of chunk
                
                for chunk_b in range(0, chromo_chunks[chr_b]+1):
                    file_name = TEMP_FILE_FORMAT.format(chr_a, chr_b, chunk_a, chunk_b)
                    chunk_path = join(temp_dir, file_name)
                    
                    if not os.path.exists(chunk_path):
                        continue
                    
                    start_b = chunk_b * CHUNK_SIZE
                    end_b = min(start_b + CHUNK_SIZE, chromo_sizes[chr_b])
                    n_b = int(math.ceil((end_b-start_b)/bin_size))
                    o_b = chromo_offsets[chr_b] + int(start_b/bin_size)
                                       
                    contact_mat = np.zeros(n_a*n_b, dtype=np.int32)
                    indices1 = np.repeat(np.arange(o_a, o_a+n_a), n_b) # Slow indices ; bin idx for chr_1
                    indices2 = np.tile(np.arange(o_b, o_b+n_b), n_a) # Fast indices ; bin id for chr_2
 
                    with io.open_file(chunk_path) as in_file_obj:
                        for i, line in enumerate(in_file_obj):

                            # Restriction fragment ENDS
                            read_id, pos_a, pos_b, strand_a, strand_b = line.split()
                            a = int((int(pos_a)-start_a)//bin_size) # Offset in chunk
                            b = int((int(pos_b)-start_b)//bin_size)
                            j = a * n_b + b
                        
                        print(f'INFO: .. chunk {chunk_a}:{chunk_b} contacts {i+1:9,d}', end='\r')
                        t += i + 1
                    
                    # Do not store zeros
                    idx = np.flatnonzero(contact_mat)
                    contact_mat = contact_mat[idx]
                    indices1 = indices1[idx]
                    indices2 = indices2[idx]
                        
                    # Global index pos
                    p2 = p1 + len(contact_mat)
 
                    hdf_pixels_bin1_id.write_direct(indices1, dest_sel=np.s_[p1:p2])
                    hdf_pixels_bin2_id.write_direct(indices2, dest_sel=np.s_[p1:p2])
                    hdf_pixels_count.write_direct(contact_mat, dest_sel=np.s_[p1:p2])
                    
                    p1 = p2

            print(f'INFO: .. {t+1:,} contacts in {time.time()-t0:.2} s       ')

    # index of chromos first bin
    hdf_indexes_chrom_offset = hdf_indexes.create_dataset('chrom_offset', (n_chromos+1), dtype='i8', compression="gzip", data=chromo_idx_bins)
    
    # index of first occurrence of bin in pixels
    hdf_indexes_bin1_offset = hdf_indexes.create_dataset('bin1_offset', (n_bins+1), dtype='i8', compression="gzip", data=bin_idx_pixels)

    hdf.flush()
    
    if i % REPORT_INTERVAL == 0:
        print(f'INFO: .. {i:,} lines read')
   
   
 
def ncc_export(ncc_paths, out_format, file_extra=None, kb_bin_sizes=[DEFAULT_BIN_SIZE]):
     
    ## only export bins with data

    join = os.path.join

    bin_sizes = sorted([int(1e3 * x) for x in kb_bin_sizes])
    
    for ncc_path in ncc_paths:
        print(f'INFO: Processing NCC file {ncc_path}')
        
        if file_extra:
            file_ext = file_extra
        else:
            file_ext = '.' + out_format
        
        if ncc_path.lower().endswith('.ncc.gz'):
            file_root = ncc_path[:-7]
        elif ncc_path.lower().endswith('.ncc'):
            file_root = ncc_path[:-4]
        else:
            file_root = ncc_path
        
        out_path = f'{file_root}{file_ext}'     
        
        if os.path.exists(out_path):
            print(f'WARNING: Previous {out_path} file will be overwritten')
            os.unlink(out_path)
        
        chromo_sizes, temp_dir = _preprocess_split_ncc(ncc_path)
        chromos = util.sort_chromosomes(chromo_sizes.keys())
        
        if out_format == 'pairs':
            
            with io.open_file(out_path, 'w') as out_file_obj:
                write = out_file_obj.write
                write('## pairs format v1.0\n')
                write('#columns: readID chr1 pos1 chr2 pos2 strand1 strand2\n')
                
                for k, chr_1 in enumerate(chromos):
                    write(f'#chromsize: {chr_1} {chromo_sizes[chr_1]}\n')
                    
                    for chr_2 in chromos[k:]:
                        if chr_1 > chr_2:
                             chr_a, chr_b = chr_2, chr_1
                        else:
                             chr_a, chr_b = chr_1, chr_2
                        
                        print(f'INFO: Writing {chr_a} {chr_b}')
                        
                        for chunk_a in range(0, int(chromo_sizes[chr_a]/CHUNK_SIZE)+1):
                            for chunk_b in range(0, int(chromo_sizes[chr_b]/CHUNK_SIZE)+1):

                               file_name = TEMP_FILE_FORMAT.format(chr_a, chr_b, chunk_a, chunk_b)
                               chunk_path = join(temp_dir, file_name)
 
                               if not os.path.exists(chunk_path):
                                   continue
                               
                               print(f'INFO: Chunk {chunk_a} {chunk_b}')
  
                               with io.open_file(chunk_path) as in_file_obj:
                                   for i, line in enumerate(in_file_obj):
 
                                       if i % REPORT_INTERVAL == 0:
                                           print(f'INFO: .. {i:,}', end='\r')

                                       # Restriction fragment ENDS
                                       read_id, pos_a,  pos_b, strand_a, strand_b = line.split()
 
                                       out_line = f'{read_id}\t{chr_a}\t{pos_a}\t{chr_b}\t{pos_b}\t{strand_a}\t{strand_b}\n'
                                       write(out_line)
 
                                   if i % REPORT_INTERVAL == 0:
                                       print(f'INFO: .. {i:,}')
 

        elif out_format == 'cool':
            hdf = h5py.File(out_path, 'w', rdcc_nbytes=8 * (1024*1024))
            print(f'INFO: Bin size {bin_sizes[0]:,}')
            _add_ncc_to_cool(ncc_path, hdf, bin_sizes[0], chromos, chromo_sizes, temp_dir)
            hdf.close()
 
        elif out_format == 'mcool':
        
            hdf = h5py.File(out_path, 'w', rdcc_nbytes=8 * (1024*1024))
            hdf.attrs['format'] = 'HDF5::MCOOL'
            hdf.attrs['format-version'] = 2
            hdf.attrs['bin-type'] = 'fixed'
                        
            hdf_resolutions = hdf.create_group('resolutions')
            
            for bin_size in bin_sizes:
                print(f'INFO: Bin size {bin_size:,}')
                hdf_group = hdf_resolutions.create_group(str(bin_size))
                _add_ncc_to_cool(ncc_path, hdf_group, bin_size, chromos, chromo_sizes, temp_dir)

            
            hdf.close()
        
        shutil.rmtree(temp_dir)
        print(f'INFO: Wrote {out_path}')
        

if __name__ == '__main__':

    from argparse import ArgumentParser
 
    epilog = 'For further help email tstevens@mrc-lmb.cam.ac.uk'

    arg_parse = ArgumentParser(prog=PROG_NAME, description=DESCRIPTION,
                               epilog=epilog, prefix_chars='-', add_help=True)

    arg_parse.add_argument(nargs='+', metavar='NCC_FILE', dest='i',
                           help='One or more input NCC format Hi-C contact Files, as generated by NucProcess. These file may be gzipped, carrying the .gz file extension.')

    arg_parse.add_argument(metavar='OUT_FORMAT', dest='f',
                           help=f'Hi-C file format in which to save data. Available: {AVAIL_FORMATS}.')

    arg_parse.add_argument('-e', metavar='FILE_EXT', default=None,
                           help='Optional suffix and file extension for naming output files based in the input files. If not specified, file extensions will be based on the file format.')
    
    mcool_default = ", ".join([f'{x:.1f}' for x in DEFAULT_BIN_SIZE_MCOOL])
    arg_parse.add_argument('-b', '--bin-size(s)', default=None, nargs='+', metavar='KB_BIN_SIZE', type=float, dest="b",
                           help=f'Binned sequence region size (the resolution) for the Hi-C contact map, in kilobases. Default for "cool" format {DEFAULT_BIN_SIZE}. For .mcool format multiple values are permitted. Default for "mcool" format {mcool_default}.')
 
    args = vars(arg_parse.parse_args())

    ncc_paths = args['i']
    out_format = args['f'].lower()
    file_extra = args['e']
    bin_sizes = args['b'] or []
    
    if out_format not in FORMATS:
        print(f'Output format (-f option) must be one of: {AVAIL_FORMATS}')
        sys.exit(0)
    
    
    if bin_sizes:
        n_resol = len(bin_sizes)
        
        if out_format == 'pairs':
            print(f'WARNING: Sequence bin size/resolution (-b option) will be ignored for "pairs" format output')
 
        elif (out_format == 'cool') and (n_resol > 1):
            print(f'FAILURE: Sequence bin size/resolution (-b option) carries multiple values, but "cool" format is specified. Maybe use "mcool" format.')
            sys.exit(0)
 
        elif (out_format == 'mcool') and (n_resol == 1):
            print(f'INFO: Sequence bin size/resolution (-b option) carries a single value, but multi-resolution "mcool" format is specified. Resolutions will be automatically set by doubling to <= 100 kb.')
            
            while bin_sizes[-1] < 100:
                bin_sizes.append(bin_sizes[-1] * 2)
            
    else: # Defaults
       if out_format == 'cool':
           bin_sizes = [DEFAULT_BIN_SIZE,]
    
       if out_format == 'mcool':
           bin_sizes = DEFAULT_BIN_SIZE_MCOOL
   
    ncc_export(ncc_paths, out_format, file_extra, bin_sizes)
