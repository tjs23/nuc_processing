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
import numpy as np

from collections import defaultdict
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

PROG_NAME = 'nuc_process_report'
VERSION = '1.0.0'
DESCRIPTION = 'Tool to present sequence processing reports from nuc_process in PDF format, given a JSON input'

def _format_table_data(d, max_nchar=70):

  data = []
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

    data.append(row)
 
  col_widths = defaultdict(int)
  row_heights = defaultdict(int)
  
  for i, row in enumerate(data):
    ht = 1
    for j, text in enumerate(row):
      text = text.strip()
      
      if len(text) > max_nchar:
        parts = []
        
        for part in text.split():
          if len(part) > max_nchar:
            while part:
              parts.append(part[:max_nchar])
              part = part[max_nchar:]
            
          else:
            parts.append(part)
        
        n = 0
        parts2 = [parts[0]]
        for part in parts[1:]:
          if len(parts2[-1]) + len(part) < max_nchar-1:
            parts2[-1] = parts2[-1] + ' ' + part
          
          else:
            parts2.append(part)
        parts = parts2
        

        text = '\n'.join(parts)
      
      nl = text.count('\n') + 1
      ht = max(nl, ht)
      wd = max_nchar if nl > 1 else len(text)
      col_widths[j] = max(col_widths[j], wd)
      row[j] = text
       
    row_heights[i] = ht
  
  char_height = sum(row_heights.values()) + 1 # Title
  char_width = sum(col_widths.values())
  
  return data, char_width, char_height, col_widths, row_heights


def _table(figure, y_start, title, table_data,  table_color, fontsize=9, x_start=0.05, dx=0.4, linespacing=1.5):
  
  width, height = figure.get_size_inches()
  
  data, char_width, char_height, col_widths, row_heights = _format_table_data(table_data)
    
  # points * inches_scale / inches_height  
  dy = linespacing * (char_height * fontsize * 0.0138889)/height
  y0 = y_start - dy
  
  ax = figure.add_axes([x_start, y0, dx, dy])
  
  y = 0.0
  t = ax.text(0.0, y, title, fontweight='bold', fontsize=fontsize, va='top')
  y += 1

  ax.hlines(y, 0, char_width, color='#B0B0B0', linewidth=1.5)
  y += 0.5
      
  for i, row in enumerate(data):
    sz = fontsize
    
    for j, text in enumerate(row):
      if j :
        x += col_widths[j-1]
      else:
        x = 0.0
      
      parts = text.split('\n')
      k = len(parts)
      t = ax.text(x, y, text, color=table_color, fontsize=fontsize, va='top', linespacing=linespacing)

      sz = max(sz, k*fontsize)

    y += row_heights[i]
    
  ax.set_xlim(0.0, char_width)
  ax.hlines(char_height+0.5, 0, char_width, color='#B0B0B0', linewidth=1.5)
  ax.set_ylim(char_height+1.0, 0.0)
  ax.set_axis_off()
  
  return ax, y0, dy
  
  
def _pie_values(data, names):

  sizes = []
  names_2 = []

  data = dict(data)

  for key in names:
    val = data[key]

    if isinstance(val, (tuple, list)):
      val, n = val
    sizes.append(val)
    name = key.replace('_', ' ')
    name = name[0].upper() + name[1:]
    names_2.append(name)

  if sizes and sum(sizes):
    return sizes, names_2

  else:
    return [100], ['No data']


def _pie_label(x):

  if x > 10.0:
    return "{:.1f}%".format(x)
  else:
    return ''

def _pie_chart(ax, stats, labels, colors):
  
  vals, labels = _pie_values(stats, labels)
  wedges, texts, autotexts = ax.pie(vals, autopct=lambda percent: _pie_label(percent),
                                    wedgeprops = {'linewidth': 0.1, 'edgecolor':'#A0A0A0'},
                                    center=(0.5, 0.5), radius=0.5, frame=True,
                                    labels=None, colors=colors, textprops=dict(color="w", fontsize=7))   
  ax.legend(wedges, labels, loc="center left", bbox_to_anchor=(1, 0, 0.5, 1), fontsize=8)
  ax.set_axis_off()

COLORS = ['#0090FF','#D0D000','#FF0000','#B0B0B0',
          '#E040E0','#FF8000','#00BBBB','#8000E0',
          '#008000','#70B070','#FFB0B0','#800000']

WATERMARK = 'nuc_process_report'

def nuc_process_report(json_stat_path, out_pdf_path=None, screen_gfx=False, fig_width=8.0, dpi=300, table_color='#006464'):
  
  with open(json_stat_path) as file_obj:
    stat_dict = json.load(file_obj)
  
  second_genome = 'map_3' in stat_dict
    
  if screen_gfx:
    pdf = None
  else:
    if not out_pdf_path:
      out_pdf_path = os.path.splitext(json_stat_path)[0] + '.pdf'
      
    pdf = PdfPages(out_pdf_path)
  
  # Clip, align, pair
   
  if len(stat_dict['command']) == 1:
    version = '1.1.0'
  else:
    version = stat_dict['command'][1][1]
  
  if second_genome:
    aspect = 1.6180339887
  else:
    aspect = 1.1

  pad = 0.025
  mid = 0.48
    
  fig = plt.figure()
  fig.set_size_inches(fig_width, fig_width * aspect)
  
  ax = fig.add_axes([0.05, 0.96, 0.9, 0.04])
  ax.text(0.0, 0.0, 'nuc_process version %s report' % version, fontweight='bold',
          fontsize=13, color=table_color)
  ax.set_axis_off()
  
  y  = 0.94
  ax, y, dy = _table(fig, y, 'Clipping reads 1', stat_dict['clip_1'], table_color)
  
  ax = fig.add_axes([mid+0.02, y, 0.4, 0.87*dy])  
  hist, hist2, edges = stat_dict['re1_pos_1']
  y_max = max(2.0, max(hist), max(hist2))
  
  ax.plot(edges, hist, color=COLORS[2], alpha=0.5, label='Unligated')
  ax.plot(edges, hist2, color=COLORS[0], alpha=0.5, label='Ligated')
  ax.set_xlabel('Read RE1 site position', fontsize=8, labelpad=2)
  ax.set_ylabel('% Reads', fontsize=8, labelpad=2)
  ax.tick_params(axis='both', which='both', labelsize=7, pad=2)
  ax.legend(fontsize=7)
  ax.set_ylim((0.0, y_max))
  ax.set_xlim((-1, max(edges)+1))
  
  y -= pad
  ax, y, dy = _table(fig, y, 'Clipping reads 2',stat_dict['clip_2'], table_color)

  ax = fig.add_axes([mid+0.02, y, 0.4, 0.87*dy])  
  hist, hist2, edges = stat_dict['re1_pos_2']
  y_max = max(2.0, max(hist), max(hist2))
  
  ax.plot(edges, hist, color=COLORS[2], alpha=0.5, label='Unligated')
  ax.plot(edges, hist2, color=COLORS[0], alpha=0.5, label='Ligated')
  ax.set_xlabel('Read RE1 site position', fontsize=8, labelpad=2)
  ax.set_ylabel('% Reads', fontsize=8)
  ax.tick_params(axis='both', which='both', labelsize=7, pad=2)
  ax.legend(fontsize=7)
  ax.set_ylim((0.0, y_max))
  ax.set_xlim((-1, max(edges)+1))

  sam_stats1 = stat_dict['map_1']
  sam_stats2 = stat_dict['map_2']
  
  sam_stats1.append(('primary_strand', stat_dict['primary_strand'][0]))
  sam_stats2.append(('primary_strand', stat_dict['primary_strand'][1]))
 
  y -= pad
  ax, y, dy = _table(fig, y,  'Genome alignment reads 1', sam_stats1, table_color)
      
  ax = fig.add_axes([mid, y, dy*aspect, dy])
  _pie_chart(ax, sam_stats1, ['unique','ambiguous','unmapped'], COLORS)
     
  y -= pad 
  ax, y, dy = _table(fig, y, 'Genome alignment reads 2', sam_stats2, table_color)
  
  ax = fig.add_axes([mid, y, dy*aspect, dy])
  _pie_chart(ax, sam_stats2, ['unique','ambiguous','unmapped'], COLORS)
  
  if second_genome:
    sam_stats3 = stat_dict['map_3']
    sam_stats4 = stat_dict['map_4']
    sam_stats3.append(('primary_strand', stat_dict['primary_strand'][2]))
    sam_stats4.append(('primary_strand', stat_dict['primary_strand'][3]))
    
    y -= pad
    ax, y, dy = _table(fig, y, 'Genome alignment 2 reads 1', sam_stats3, table_color)
   
    ax = fig.add_axes([mid, y, dy*aspect, dy])
    _pie_chart(ax, sam_stats3, ['unique','ambiguous','unmapped'], COLORS)

    y -= pad
    ax, y, dy = _table(fig, y, 'Genome alignment 2 reads 2', sam_stats4, table_color)

    ax = fig.add_axes([mid, y, dy*aspect, dy])
    _pie_chart(ax, sam_stats4, ['unique','ambiguous','unmapped'], COLORS)
  
  if second_genome:
    y -= pad
    pair_stats = stat_dict['pair']
    yax, y, dy = _table(fig, y, 'Pairing reads', pair_stats, table_color)
 
    ax = fig.add_axes([mid, y, dy*aspect, dy])
    _pie_chart(ax, pair_stats, ['unique','position_ambiguous','unmapped_end','genome_ambiguous'], COLORS)
  else:
    y -= pad
    pair_stats = stat_dict['pair']
    ax, y, dy = _table(fig, y, 'Pairing reads', pair_stats, table_color)
 
    ax = fig.add_axes([mid, y, dy*aspect, dy])
    _pie_chart(ax, pair_stats, ['unique','ambiguous','unmapped_end'], COLORS)
  
  ax.text(0.01, 0.01, WATERMARK, color='#B0B0B0', fontsize=8, transform=fig.transFigure)
   
  if pdf:
    pdf.savefig(dpi=dpi)
  else:
    plt.show()  
   
  # Filter
  

  redundancy_stats = stat_dict.get('dup')
  promiscuity_stats = stat_dict.get('promsic')
  
  if promiscuity_stats:
    aspect = 1.61803398875
  else:
    aspect = 1.25
    
  pad = 0.01
  
  fig = plt.figure()
  fig.set_size_inches(fig_width, fig_width * aspect)
  
  y = 0.94
 
  filter_stats = stat_dict['filter']
  
  for x in filter_stats:
    if x[0] == 'excluded_ambig_group':
      x[0] = 'excluded_group'
      
  ax, y, dy = _table(fig, y, 'Filtering pairs', filter_stats, table_color)
  
  ax = fig.add_axes([mid, y+0.1*dy , 0.75*dy*aspect, 0.75*dy])
  
  labels = ['accepted','internal_re1','adjacent_re1','circular_re1',
            'overhang_re1', 'too_close','too_small','too_big',
            'internal_re2','no_end_re2', 'unknown_contig', 'excluded_group']

  _pie_chart(ax, filter_stats, labels, COLORS)
  
  dy = 0.2
  y -= dy
  y -= pad
  y -= pad
  
  hist1, hist2, hist3, edges13, hist4, hist5, edges45 = [np.array(x) for x in stat_dict['frag_sizes']]
  lim = hist1.max()/2e2
  idx = (hist1 > lim).nonzero()[0]
  i = idx[1]-1
  j = idx[-1]+1
  
  ax = fig.add_axes([0.08, y, 0.4, dy])
  
  edges13 = 10.0**edges13[i:j]
  ax.semilogx(edges13, hist1[i:j], color=COLORS[3], alpha=0.5, label='All')
  ax.semilogx(edges13, hist2[i:j], color=COLORS[2], alpha=0.5, label='$Cis$ accepted')
  ax.semilogx(edges13, hist3[i:j], color=COLORS[0], alpha=0.5, label='$Trans$ accepted')
  ax.legend(fontsize=8)
  ax.set_xlabel('Fragment size $log_{10}$(bp)', fontsize=8, labelpad=2)
  ax.set_ylabel('% Pairs', fontsize=8, labelpad=2)
  ax.tick_params(axis='both', which='both', labelsize=7, pad=2)
  
  ax = fig.add_axes([mid+0.05, y, 0.4, 0.20])

  ax.plot(edges45, hist4, color='#FF0000', alpha=0.5, label='$+$ strand')
  ax.plot(edges45, hist5, color='#0090FF', alpha=0.5, label='$-$ strand')
  ax.legend(fontsize=8)
  ax.set_xlabel('RE1 site separation (bp)', fontsize=8, labelpad=2)
  ax.set_ylabel('Reads/bp', fontsize=8, labelpad=2)
  ax.tick_params(axis='both', which='both', labelsize=7, pad=2)
  y -= 0.04
  
  if redundancy_stats is not None:
    y -= pad
    ax, y, dy = _table(fig, y, 'Duplicate removal', redundancy_stats, table_color)

    ax = fig.add_axes([mid, y, dy*aspect, dy])
    _pie_chart(ax, redundancy_stats, ['redundant','unique'], COLORS)

  
  if promiscuity_stats is not None:
    y -= pad
    ax, y, dy = _table(fig, y, 'Promiscuous pair removal', promiscuity_stats, table_color)

    ax = fig.add_axes([mid, y, dy*aspect, dy])
    _pie_chart(ax, promiscuity_stats, ['clean','resolved','promiscuous'], COLORS)
  
  
  y -= pad
  final_stats = stat_dict['final']
  ax, y, dy = _table(fig, y,  'Final output', final_stats, table_color)
  
  ax = fig.add_axes([mid, y, dy*aspect, dy])
  
  names = ['trans','cis_far','cis_near']
  
  if second_genome:
    names.append('homolog_trans')
  
  _pie_chart(ax, final_stats, names, COLORS)
 
  ax.text(0.01, 0.01, WATERMARK, color='#B0B0B0', fontsize=8, transform=fig.transFigure)
  
  if pdf:
    pdf.savefig(dpi=dpi)
  else:
    plt.show()  
 
  
  # Inputs
 
  input_stats = stat_dict['general']
  if len(stat_dict['command']) > 1:
    input_stats.append(stat_dict['command'][1])
  input_stats.append(stat_dict['command'][0])
  
  null, char_width, char_height, col_widths, row_heights = _format_table_data(input_stats)
  
  fig = plt.figure()
  fig.set_size_inches(fig_width, 1.5 * (3.0+char_height) * 9.0 * 0.0138889)
  
  ax, y, dy = _table(fig, 0.94, 'Input parameters', input_stats, table_color, dx=0.9)

  ax.text(0.01, 0.01, WATERMARK, color='#B0B0B0', fontsize=8, transform=fig.transFigure)
   
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

  arg_parse.add_argument('-o', '--pdf-out', metavar='OUT_FILE', nargs='+', default=None, dest='o',
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
