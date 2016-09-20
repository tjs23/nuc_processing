import numpy as np
import os, sys

from math import ceil
from NucSvg import SvgDocument

PROG_NAME = 'nuc_contact_map'
VERSION = '1.0.0'
DESCRIPTION = 'Chromatin contact (NCC format) Hi-C contact map display module for Nuc3D and NucTools'

# TBC:
# - ambiguity filter

def info(msg, prefix='INFO'):

  print('%8s : %s' % (prefix, msg))


def warn(msg, prefix='WARNING'):

  print('%8s : %s' % (prefix, msg))


def fatal(msg, prefix='%s FAILURE' % PROG_NAME):

  print('%8s : %s' % (prefix, msg))
  sys.exit(0)


def _color_func_white(matrix, color):
  
  n, m = matrix.shape
  color_matrix = np.zeros((n, m, 3), float)
  color = np.array(color)
  rgb_0 = np.array([255.0, 255.0, 255.0])
  rgb_b = np.array([0.0, 0.0, 0.0])  
  
  for i in range(n):
    for j in range(m):
      f = matrix[i,j]
      
      if f > 0:
        g = 1.0 - f
        color_matrix[i,j] = (g * color) + (f * rgb_b)
         
      else:
        color_matrix[i,j] = rgb_0

  color_matrix = np.clip(color_matrix, 0, 255)
  color_matrix = np.array(color_matrix, dtype=np.uint8)

  return color_matrix


def _color_func_black(matrix, color):
  
  n, m = matrix.shape
  color_matrix = np.zeros((n, m, 3), float)
  color = np.array(color)
  rgb_0 = np.array([0.0, 0.0, 0.0])
  rgb_b = np.array([255.0, 255.0, 255.0])
  
  for i in range(n):
    for j in range(m):
      f = matrix[i,j]
      
      if f > 0:
        g = 1.0 - f
        color_matrix[i,j] = (g * color) + (f * rgb_b)
         
      else:
        color_matrix[i,j] = rgb_0

  color_matrix = np.clip(color_matrix, 0, 255)
  color_matrix = np.array(color_matrix, dtype=np.uint8)

  return color_matrix


def nuc_contact_map(ncc_path, svg_path=None, svg_width=700, bin_size=5, black_bg=False, color=None, font=None, font_size=12, line_width=1):
    
  bin_size = int(bin_size * 1e6)
  chromo_limits = {}  
  
  # Load NCC data
  contacts = {}
  with open(ncc_path) as in_file_obj:
    for line in in_file_obj:
      chr_a, f_start_a, f_end_a, start_a, end_a, strand_a, chr_b, f_start_b, f_end_b, start_b, end_b, strand_b, ambig_group, pair_id, swap_pair = line.split()
      
      f_start_a = int(f_start_a)
      f_end_a = int(f_end_a)
      start_a = int(start_a)
      end_a = int(end_a)
      f_start_b = int(f_start_b)
      f_end_b = int(f_end_b)
      start_b = int(start_b)
      end_b  = int(end_b)
      
      if chr_a in chromo_limits:
        s, e = chromo_limits[chr_a]
        chromo_limits[chr_a] = (min(s, f_start_a), max(e, f_end_a))
      else:
        chromo_limits[chr_a] = (f_start_a, f_end_a)

      if chr_b in chromo_limits:
        s, e = chromo_limits[chr_b]
        chromo_limits[chr_b] = (min(s, f_start_b), max(e, f_end_b))
      else:
        chromo_limits[chr_b] = (f_start_b, f_end_b)
      
      # Store data in sorted chromo order
      if chr_a > chr_b:
        chr_a, chr_b = chr_b, chr_a
        start_a, start_b = start_b, start_a
        end_a, end_b = end_b, end_a
      
      chr_pair = (chr_a, chr_b)
      if chr_pair in contacts:
        contact_list = contacts[chr_pair]
      
      else:
        contact_list = []
        contacts[chr_pair] = contact_list
        
      contact_list.append((start_a, end_a, start_b, end_b))

  if not chromo_limits:
    fatal('No chromosome contact data read')
  
  # Get sorted chromosomes    
  chromos = []
  for chromo in chromo_limits:
    if chromo.upper().startswith('CHR'):
      c = chromo[3:]
    else:
      c = chromo
    
    try:
      key = '%09d' % int(c)
    except ValueError as err:
      key = c
 
    chromos.append((key, chromo))
   
  chromos.sort()
  chromos = [x[1] for x in chromos]
  
  # Get chromosome matrix index ranges
  grid = []
  chromo_offsets = {}
  chromo_spans = {}
  n = 0 
  for chromo in chromos: # In display order
    s, e = chromo_limits[chromo]
    a = int(s/bin_size)
    b = int(ceil(e/float(bin_size))) + 1
    span = b-a
    chromo_offsets[chromo] = s, n # Start bp, start bin index
    chromo_spans[chromo] = span
    n += span
    grid.append(n)# At chromosome edge
  
  grid.pop() # Don't need last edge
  # Fill contact map matrix
  data = np.zeros((n, n), float)
  
  if svg_path:
    info('Contact map size %d x %d' % (n, n))
  
  for i, chr_1 in enumerate(chromos):
    for chr_2 in chromos[i:]:
      
      if chr_1 > chr_2:
        chr_a, chr_b = chr_2, chr_1
      else:
        chr_a, chr_b = chr_1, chr_2
        
      contact_list = contacts.get((chr_a, chr_b))
      
      if contact_list is None: # Nothing for this pair: common for single-cell Hi-C
        continue
      
      s_a, n_a = chromo_offsets[chr_a]
      s_b, n_b = chromo_offsets[chr_b]
      
      for start_a, end_a, start_b, end_b in contact_list:
        a1 = n_a + int((start_a-s_a)/bin_size)
        a2 = n_a + int((end_a-s_a)/bin_size)
        b1 = n_b + int((start_b-s_b)/bin_size)
        b2 = n_b + int((end_b-s_b)/bin_size)
        
        fa = 1.0 / (a2 - a1 + 1.0)
        fb = 1.0 / (b2 - b1 + 1.0)
        
        for a in range(a1, a2+1):
          for b in range(b1, b2+1):
            f =  fa * fb
            data[a, b] += f
            data[b, a] += f
  
  data = np.log(data+1.0)
  
  chromo_labels = []
  for chromo in chromos:
    pos = chromo_offsets[chromo][1] + chromo_spans[chromo]/2
    
    if chromo.upper().startswith('CHR'):
      chromo = chromo[3:]
    
    chromo_labels.append((pos, chromo))

  
  # Make SVG  
  offset = int(0.1 * n)  
  
  svg_doc = SvgDocument()  
  
  
  if black_bg:
    if color:
      color = svg_doc._hex_to_rgb(color)
    else:
      color = [55.0, 55.0, 110.0]
    
    color_func = lambda x, c=color: _color_func_black(x, c)
    grid_color = '#303030'
    
  else:
    if color:
      color = svg_doc._hex_to_rgb(color)
    else:
      color = [200.0, 200.0, 255.0]
   
    color_func = lambda x, c=color: _color_func_white(x, c)
    grid_color = '#C0C0C0'
  
  w = svg_width - 2 * offset
  svg_doc.density_matrix(data, w, w, x_grid=grid, y_grid=grid,
                         x_labels=chromo_labels, y_labels=chromo_labels,
                         x_axis_label='Chromosome', y_axis_label='Chromosome',
                         grid_color=grid_color, font=font,
                         font_size=font_size, line_width=line_width,
                         plot_offset=(offset, offset), color_func=color_func,
                         value_range=None, scale_func=None)
  
  svg_doc.text(os.path.basename(ncc_path), (offset, offset/2))

  if svg_path == '-':
    print svg_doc.svg()
  
  elif svg_path:
    svg_doc.write_file(svg_path, svg_width, svg_width)

  else:
    return svg_doc.svg()


if __name__ == '__main__':
 
  from argparse import ArgumentParser
  
  epilog = 'For further help email tjs23@cam.ac.uk or wb104@cam.ac.uk'
  
  arg_parse = ArgumentParser(prog=PROG_NAME, description=DESCRIPTION,
                             epilog=epilog, prefix_chars='-', add_help=True)

  arg_parse.add_argument('-i', metavar='NCC_FILE',
                         help='Input NCC format chromatin contact file')

  arg_parse.add_argument('-o', metavar='SVG_FILE',
                         help='Optional output name for SVG format contact map output')

  arg_parse.add_argument('-w', default=700, metavar='SVG_WIDTH',
                         type=int, help='SVG document width')

  arg_parse.add_argument('-s', default=5, metavar='BIN_SIZE',
                         help='Sequence region size represented by each small square (the resolution) in megabases. Default is 5 kb')

  arg_parse.add_argument('-b', default=False, action='store_true',
                         help='Specifies that the contact map should have a black background (default is white)')

  arg_parse.add_argument('-c', nargs=1, metavar='RGB_COLOR',
                         help='Optional main color for the contact points as a 24-bit hexidecimal RBG code e.g. "#0080FF" (with quotes)')

  args = vars(arg_parse.parse_args())
  
  ncc_path  = args['i']
  svg_path  = args['o']
  svg_width = args['w']
  bin_size  = args['s']
  black_bg  = args['b']
  color     = args['c']
  
  if not ncc_path:
    import sys
    arg_parse.print_help()
    sys.exit(1)
    
  if color and isinstance(color, list):
    color = color[0]
  
  if not svg_path:
    svg_path = os.path.splitext(ncc_path)[0] + '_contact_map.svg'
    
  nuc_contact_map(ncc_path, svg_path, svg_width, bin_size, black_bg, color)

