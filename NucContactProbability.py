import sys
import os
import numpy as np
import NucSvg

from matplotlib import pyplot as plt


PROG_NAME = 'nuc_contact_probability'
VERSION = '1.0.0'
DESCRIPTION = 'Chromatin contact (NCC format) probability vs sequence separation module for Nuc3D and NucTools'


def info(msg, prefix='INFO'):

  print('%8s : %s' % (prefix, msg))


def warn(msg, prefix='WARNING'):

  print('%8s : %s' % (prefix, msg))


def fatal(msg, prefix='%s FAILURE' % PROG_NAME):

  print('%8s : %s' % (prefix, msg))
  sys.exit(0)
  
  
def plot_contact_probability_seq_sep(ncc_paths, svg_path, bin_size=100, svg_width=800):
 
  bin_size *= 1e3
  n_files = len(ncc_paths)
 
  chromo_limits = {} 
  contacts = {}
  seq_seps_all = []
  seq_seps = {}
  weights_all = []
  weights = {}
   
  for ncc_path in ncc_paths:
    seq_seps[ncc_path] = []
    weights[ncc_path] = []
  
    # Load cis NCC data
    with open(ncc_path) as in_file_obj:
      for line in in_file_obj:
        chr_a, f_start_a, f_end_a, start_a, end_a, strand_a, chr_b, f_start_b, f_end_b, start_b, end_b, strand_b, ambig_group, pair_id, swap_pair = line.split()
        
        if chr_a != chr_b:
          continue
         
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
          chromo_limits[chr_a] = (min(s, f_start_a, f_start_b), max(e, f_end_a, f_end_b))
        else:
          chromo_limits[chr_a] = (min(f_start_a, f_start_b), max(f_end_a, f_end_b))
          
        if chr_a in contacts:
          contact_list = contacts[chr_a]
 
        else:
          contact_list = []
          contacts[chr_a] = contact_list
 
        if strand_a > 0:
          p_a = f_end_a
        else:
          p_a = f_start_a
 
        if strand_b > 0:
          p_b = f_end_b
        else:
          p_b = f_start_b
 
        contact_list.append((p_a, p_b))

    for chromo in contacts:
      
      contact_array = np.array(contacts[chromo], np.int32)
      
      seps = abs(contact_array[:,0]-contact_array[:,1])

      indices = seps.nonzero()
      seps    = seps[indices]
      
      p_start, p_end = chromo_limits[chromo]
      size = float(p_end-p_start+1)
      
      prob = (size/(size-seps)).tolist() # From fraction of chromosome that could give rise to each separation
      seps = seps.tolist()
            
      seq_seps[ncc_path] += seps
      seq_seps_all += seps
      
      weights[ncc_path] += prob
      weights_all += prob

  seq_seps_all = np.array(seq_seps_all)
       
   
  
  log_max  = (int(2 * np.log10(seq_seps_all.max()))/2.0)
  log_min  = np.log10(bin_size)
  x_limit = 10.0 ** log_max 
  
  num_bins = (x_limit-bin_size)/bin_size
  bins = np.linspace(bin_size, x_limit, num_bins)
      
  opacities = []
  data_lists = []
  colors = []
  names = []
  
  f, ax = plt.subplots()
  
  ax.set_title('Contact sequence separations')
   
  ax.set_alpha(0.5)  
  
  if n_files > 1:
    for i, ncc_path in enumerate(seq_seps.keys()):
      data = np.array(seq_seps[ncc_path])

      hist, edges = np.histogram(data, bins=bins, weights=weights[ncc_path], normed=True)
      
      idx = hist.nonzero()
 
      hist = hist[idx]
      edges = edges[idx]
 
      x_data = np.log10(edges)
      y_data = np.log10(hist)
      
      if i == 0:
        label = 'Individual datasets'
      else:
        label = None
      
      data_lists.append(zip(x_data, y_data))
      colors.append('#FF0000')
      names.append(label)
      opacities.append(0.25)
        
      ax.plot(x_data, y_data, label=label, color='#FF0000', alpha=0.25)
  
  hist, edges = np.histogram(seq_seps_all, bins=bins, weights=weights_all, normed=True) 
  idx = hist.nonzero()
  
  hist = hist[idx]
  edges = edges[idx]
  
  x_data = np.log10(edges)
  y_data = np.log10(hist) 
  
  y_min = int(2.0 * y_data.min())/2.0
  y_max = int(1.0 + 2.0 * y_data.max())/2.0
  
  if n_files > 1:
    label = 'Combined datasets'
  
  else:
    label = None
  
  data_lists.append(zip(x_data, y_data))
  colors.append('#000000')
  names.append(label)
  opacities.append(1.0)
    
  ax.plot(x_data, y_data, label=label, color='#000000', linewidth=1, alpha=1.0)
  
  x1 = x_data[2]
  y1 = y_data[2]
  dy = 0.25
  
  ax.plot([x1, x1+1.5], [y1+dy, y1+dy-1.50], color='#808080', linestyle='--', alpha=0.5)   
  ax.text(x1+1.5, y1+dy-1.50, r'$\lambda$=1.0', color='#808080', fontsize=16, alpha=0.5, va='center')

  ax.plot([x1, x1+1.5], [y1-dy, y1-dy-2.25], color='#808080', linestyle='--', alpha=0.5) 
  ax.text(x1+1.5, y1-dy-2.25, r'$\lambda$=1.5', color='#808080', fontsize=16, alpha=0.5, va='center')

  data_lists.append([(x1, y1+dy), (x1+1.5, y1+dy-1.50)])
  colors.append('#808080')
  names.append(None)
  opacities.append(0.5)
  
  data_lists.append([(x1, y1-dy), (x1+1.5, y1-dy-2.25)])
  colors.append('#808080')
  names.append(None)
  opacities.append(0.5)

  ax.legend()
  ax.set_xlabel('Sequence separation (bp)')
  ax.set_ylabel('Contact probability (%d kb bins)' % (bin_size/1000))
 
  x_range = np.arange(np.log10(bin_size), np.log10(x_limit), 0.5)
  
  ax.xaxis.set_ticks(x_range)
  ax.set_xticklabels(['$10^{%.1f}$' % x for x in x_range], fontsize=12)

  y_range = np.arange(y_min, y_max, 0.5)
  ax.yaxis.set_ticks(y_range)
  ax.set_yticklabels(['$10^{%.1f}$' % x for x in y_range], fontsize=12)
  
  ax.set_xlim((x_range[0], x_range[-1]+0.25))
  ax.set_ylim((y_min, y_max))
  
  plt.show()

  svg_doc = NucSvg.SvgDocument()
  
  pad = 100
  width   = svg_width - pad
  height  = svg_width - pad
  x_label = 'Sequence separation (bp)'
  y_label = 'Contact probability (%d kb bins)' % (bin_size/1000)
  title   = 'Contact sequence separations'
  texts = [(u'\u03bb=1.0', '#808080', (x1+1.5, y1+dy-1.50)),
           (u'\u03bb=1.5', '#808080', (x1+1.5, y1-dy-2.25))]
   
  svg_doc.graph(0, 0, width, height, data_lists, x_label, y_label,
            names=names, colors=colors,  graph_type=NucSvg.LINE_TYPE,
            symbols=None, line_widths=None, symbol_sizes=None,
            legend=(6.4, -6.0), 
            title=title, x_labels=None, plot_offset=(pad, pad),
            axis_color='black', bg_color='#F0F0F0', font=None, font_size=16, line_width=1,
            x_ticks=True, y_ticks=True, x_grid=False, y_grid=False,
            texts=texts, opacities=opacities,
            x_log_base='10', y_log_base='10')
    
  if svg_path == '-':
    print(svg_doc.svg(svg_width, svg_width))
  
  elif svg_path:
    svg_doc.write_file(svg_path, svg_width, svg_width)

  else:
    return svg_doc.svg(svg_width, svg_width)
  


if __name__ == '__main__':
 
  from argparse import ArgumentParser
  
  epilog = 'For further help email tjs23@cam.ac.uk or wb104@cam.ac.uk'
  
  arg_parse = ArgumentParser(prog=PROG_NAME, description=DESCRIPTION,
                             epilog=epilog, prefix_chars='-', add_help=True)

  arg_parse.add_argument('ncc_paths', metavar='NCC_FILE', nargs='+',
                         help='Input NCC format chromatin contact file(s). Wildcards accepted')

  arg_parse.add_argument('-o', metavar='NCC_FILE',
                         help='Output SVF format file')

  arg_parse.add_argument('-w', default=800, metavar='SVG_WIDTH', type=int,
                         help='SVG document width')

  arg_parse.add_argument('-s', default=100, metavar='KB_BIN_SIZE', type=int,
                         help='Sequence region size in kilobases for calculation of contact probabilities. Default is 100 (kb)')

  args = vars(arg_parse.parse_args())
  
  ncc_paths = args['ncc_paths']
  svg_path  = args['o']
  svg_width = args['w']
  bin_size  = args['s']
  
  
  if not ncc_paths:
    import sys
    arg_parse.print_help()
    fatal('No input NCC format files specified') 

  for ncc_path in ncc_paths:
    if not os.path.exists(ncc_path):
      fatal('NCC file could not be found at "{}"'.format(ncc_path))  
  
  if not svg_path:
    arg_parse.print_help()
    fatal('Output SVG file path not specified') 
  
  dir_name, file_name = os.path.split(svg_path)
  
  if dir_name and not os.path.exists(dir_name):
    fatal('Directory "%s" for output could not be found') 
  
  if not file_name.lower().endswith('.svg'):
    svg_path = svg_path + '.svg'
  
  plot_contact_probability_seq_sep(ncc_paths, svg_path, bin_size, svg_width)
