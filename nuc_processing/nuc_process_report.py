"""
---- COPYRIGHT ----------------------------------------------------------------

Copyright (C) 20016-2019
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

import os, json, sys, math

from collections import defaultdict
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

PROG_NAME = 'nuc_process_report'
VERSION = '1.0.0'
DESCRIPTION = 'Tool to present sequence processing reports from nuc_process in PDF format, given a JSON input'

def _format_list(d):

  l = []
  for k, v in d:
    c1 = k.replace('_', ' ')
    c1 = c1[0].upper() + c1[1:]
    
    row = [c1]
    
    if isinstance(v, (tuple, list)):
      v, n = v
      percent = 100.0 * v/float(n or 1)
      row += ['{:,}'.format(v), '{:.2f}%'.format(percent)]
      
    elif isinstance(v, int):
      row += ['{:,}'.format(v)]

    elif isinstance(v, float):
      row += ['{:.3f}'.format(v)]

    else:
      row += [v]

    l.append(row)

  return l


def _table(ax, title, data,  table_color, fontsize=9):
  
  n_rows = 1 # Title
  col_sz = defaultdict(int)
  row_sz = defaultdict(int) 
  
  if len(data[0]) > 2:
    max_nchar = 30
    
  elif len(data[0]) > 1:  
    max_nchar = 70  
 
  else:
    max_nchar = 100
  
  for i, row in enumerate(data):
    ht = 1
    for j, text in enumerate(row):
      text = text.strip()
      
      if len(text) > max_nchar:
        parts = text.split()
        
        if parts:
          n = 0
          parts2 = [parts[0]]
          for part in parts[1:]:
            if len(parts2[-1]) + len(part) < max_nchar-1:
              parts2[-1] = parts2[-1] + ' ' + part
            
            else:
              parts2.append(part)
          parts = parts2
          
        else:
          parts = []
          while text:
            parts.append(text[:max_nchar])
            text = text[max_nchar:]
      
        text = '\n'.join(parts)
      
      nl = text.count('\n') + 1
      ht = max(nl, ht)
      wd = max_nchar if nl > 1 else len(text)
      col_sz[j] = max(col_sz[j], wd)
      row[j] = text
       
    row_sz[i] = ht
    n_rows += ht
  
  n_char = sum(col_sz.values())
  
  y = 0.0
  t = ax.text(0.0, y, title, fontweight='bold', fontsize=fontsize, va='top')
  y += fontsize + 1

  ax.hlines(y, 0, n_char, color='#B0B0B0', linewidth=1.5)
  y += 3
      
  for i, row in enumerate(data):
    sz = fontsize
    
    for j, text in enumerate(row):
      if j :
        x += col_sz[j-1]
      else:
        x = 0.0
      
      k = text.count('\n') + 1
      t = ax.text(x, y, text, color=table_color, fontsize=fontsize, va='top')

      sz = max(sz, k*fontsize)

    y += sz # row_sz[i]
    
  ax.set_xlim(0.0, n_char)
  ax.hlines(y, 0, n_char, color='#B0B0B0', linewidth=1.5)
  ax.set_ylim(y+1, 0.0)
  ax.set_axis_off()
  

def nuc_process_report(json_stat_path, out_pdf_path=None, screen_gfx=False, fig_width=8.0, dpi=300, table_color='#006464'):
 
  with open(json_stat_path) as file_obj:
    stat_dict = json.load(file_obj)
    
  if screen_gfx:
    pdf = None
  else:
    if not out_pdf_path:
      out_pdf_path = os.path.splitext(json_stat_path)[0] + '.pdf'
      
    pdf = PdfPages(out_pdf_path)
  
  # Inputs
  
  fig = plt.figure()
  fig.set_size_inches(fig_width, fig_width * 1.61803398875)
  
  
  
  if len(stat_dict['command']) == 1:
    command = stat_dict['command'][0][1]
    version = '1.1.0'
  else:
    command = stat_dict['command'][0][1]
    version = stat_dict['command'][1][1]
    
  ax = fig.add_axes([0.05, 0.96, 0.9, 0.04])
  ax.text(0.0, 0.0, 'nuc_process version %s report' % version, fontweight='bold', fontsize=13, color=table_color)
  ax.set_axis_off()
 
  ax = fig.add_axes([0.05, 0.65, 0.9, 0.30])
  _table(ax, 'Input parameters', _format_list(stat_dict['general']), table_color)
  
  ax = fig.add_axes([0.05, 0.30, 0.9, 0.30])
  _table(ax, 'Command used', [[command]], table_color)
  
  
  # Clip, align, pair
  
  
  # Filter
  

  if pdf:
    pdf.savefig(dpi=dpi)
    pdf.close()
    print('Info: Written {}'.format(out_pdf_path))
  
  else:
    plt.show()  
    print('Info: Done')

  plt.close()

  
def main(argv=None):

  from argparse import ArgumentParser  
  
  if argv is None:
    argv = sys.argv[1:]

  epilog = 'For further help email tjs23@cam.ac.uk or wb104@cam.ac.uk'

  arg_parse = ArgumentParser(prog=PROG_NAME, description=DESCRIPTION,
                             epilog=epilog, prefix_chars='-', add_help=True)

  arg_parse.add_argument(metavar='JSON_FILE', nargs='+', dest='i',
                         help='One or more JSN statistics files as output from nuc_process. Wildcards accepted.')

  arg_parse.add_argument('-o', '--pdf-out', metavar='OUT_FILE', default=None, dest='o',
                         help='One or more ptional output file PDF names; to match corresponding input files. ' \
                              'Defaults based on the input JSON files will be used if not specified')

  arg_parse.add_argument('-g', '--screen-gfx', default=False, action='store_true', dest='g',
                         help='Display graphics on-screen using matplotlib, where possible and do not automatically save output.')
  
  args = vars(arg_parse.parse_args(argv))
                              
  in_paths = args['i']
  out_paths = args['o'] or []
  screen_gfx = args['g']
  
  while len(out_paths) < len(in_paths):
    out_paths.append(None)
  
  in_paths = in_paths[:len(out_paths)]
  
  for json_stat_path, out_pdf_path in zip(in_paths, out_paths):
    nuc_process_report(json_stat_path, out_pdf_path, screen_gfx)
 
  
if __name__ == "__main__":
  sys.path.append(os.path.dirname(os.path.dirname(__file__)))
  main()
