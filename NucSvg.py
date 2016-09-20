from collections import defaultdict
from math import sin, cos, pi
import colorsys, base64, cStringIO
import numpy as np
from scipy import misc

TAU = 2*pi
LINE_TYPE = 'line'
SCATTER_TYPE = 'scatter'
BAR_TYPE = 'histogram'
DEFAULT_COLORS = ['#800000','#000080',
                  '#008000','#808000',
                  '#800080','#008080',
                  '#808080','#000000',
                  '#804000','#004080']

class SvgDocument(object):

  def __init__(self, font='Arial'):

    self.font = font
    self._svg_lines = []
    
    
  def clear(self):
    
    self._svg_lines = []
  
  
  def svg(self, width, height):
    
    return self._svg_head(width, height) + self._svg_lines + self._svg_tail()
    
  
  def write_file(self, file_name, width, height):
    
    file_obj = open(file_name, 'w')
    
    file_obj.writelines(self._svg_head(width, height))
    file_obj.writelines(self._svg_lines)
    file_obj.writelines(self._svg_tail())   
    file_obj.close()
  
  
  def _svg_head(self, width, height):

    head1 = '<?xml version="1.0"?>\n'
    head2 = '<svg height="%d" width="%d" image-rendering="optimizeSpeed" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" version="1.2" baseProfile="tiny">\n' % (height,width)
    
    return head1, head2
    

  def _svg_tail(self):

    return ['</svg>\n']

 
  def group_start(self, opacity=1.0, stroke='black', width=1):
 
    data = (opacity, stroke, width)
 
    self._svg_lines += ['  <g style="fill-opacity:%.1f; stroke:%s; stroke-width:%d;">\n' % data]
    
  
  def group_end(self):
 
    self._svg_lines.append('  </g>\n')


  def poly(self, coords, color='black', fill='blue', width=1):
  
    path = ['M']
    for x, y in coords:
      path += ['%d' % x, '%d' % y, 'L']
    path[-1] = 'Z'
      
    line = '    <path d="%s" fill="%s" stroke="%s" stroke-width="%d" />\n' % (' '.join(path), fill, color, width)
    self._svg_lines += line
    

  def lines(self, coords, color='black', width=1):
  
    path = ['M']
    for x, y in coords:
      path += ['%d' % x, '%d' % y, 'L']
    del path[-1]
      
    line = '     <path d="%s" fill-opacity="0.0" stroke="%s" stroke-width="%d" />\n' % (' '.join(path), color, width)
    
    self._svg_lines += line


  def line(self, coords, color='black', line_width=1):
  
    x1, y1, x2, y2 = coords
    line = '     <line x1="%d" y1="%d" x2="%d" y2="%d" stroke="%s" stroke-width="%d" />\n' % (x1, y1, x2, y2, color, line_width)
    self._svg_lines += line
    
    
  def circle(self, center, radius, color='black', fill='#808080'):

    cx, cy = center
    x = cx + radius
    y = cy + radius
    
    line = '     <circle cx="%d" cy="%d" r="%d" stroke="%s" fill="%s" />\n' % (x,y,radius,color,fill)
    self._svg_lines += line


  def rect(self, coords, color='black', fill='#808080'):

    x0, y0, x1, y1 = coords
    
    x = min((x0, x1))
    y = min((y0, y1))
    w = abs(x1 - x0)
    h = abs(y1 - y0)
    line = '     <rect x="%d" y="%d" height="%d" width="%d" stroke="%s" fill="%s" />\n' % (x,y,h,w,color,fill)
    self._svg_lines += line


  def square(self, center, radius, color='black', fill='#808080'):

    cx, cy = center
    x = cx - radius
    y = cy - radius
    r = 2*radius
    line = '     <rect x="%d" y="%d" height="%d" width="%d" stroke="%s" fill="%s" />\n' % (x,y,r,r,color,fill)
    self._svg_lines += line


  def text(self, text, coords, anchor='start', size=16, bold=False, font=None, color=None, angle=None, vert_align=None):
    
    if not font:
      font = self.font
    
    x, y = coords
    
    attrs = ''
    
    if vert_align:
      attrs += ' dominant-baseline="%s"' % vert_align
    
    if color:
      attrs += ' fill="%s"' % color
    
    if angle:
      attrs += ' transform="rotate(%d %d %d)"' % (angle, x,y)
    
    if bold:
      attrs += ' font-weight="bold"'
    
    line = '     <text x="%d" y="%d" text-anchor="%s" font-family="%s" font-size="%d"%s>%s</text>\n' % (x,y,anchor,font,size,attrs,text)
     
    self._svg_lines += line

  def segment(self, x0, y0, r, a1, a2, color='black', fill='grey', line_width=1):
    
    a1 = a1 % TAU
    a2 = a2 % TAU
    
    if a2 < a1:
      a1, a2 = a1, a2
    
    da = (a2-a1) % TAU
    la = 0 if da <= pi else 1
    clockwise = 1
    x1 = x0 + r * cos(a1)
    y1 = y0 + r * sin(a1)
    x2 = x0 + r * cos(a2)
    y2 = y0 + r * sin(a2)
  
    path = 'M %d %d L %d %d A %d %d %d %d %d %d %d Z' % (x0, y0, x1, y1, r, r, 0, la, clockwise, x2, y2)

    line = '    <path d="%s" fill="%s" stroke="%s" stroke-width="%d" />\n' % (path, fill, color, line_width)
    self._svg_lines += line
  
  
  def image(self, x, y, w, h, data):
  
    io = cStringIO.StringIO()
    img = misc.toimage(data)
    img.save(io, format="PNG")
    base_64_data = base64.b64encode(io.getvalue())
     
    line = '     <image x="%d" y="%d" width="%d" height="%d" xlink:href="data:image/png;base64,%s" />\n' % (x, y, w, h, base_64_data)
     
    self._svg_lines += line
      
      
  def _graph_get_draw_coords(self, data_points, data_region, plot_region):

    dr0, dr1, dr2, dr3 = data_region
    pr0, pr1, pr2, pr3 = plot_region
    
    delta_x_plot = pr2 - pr0
    delta_y_plot = pr3 - pr1
    delta_x_data = dr2 - dr0
    delta_y_data = dr3 - dr1
  
    ppv_x = delta_x_plot/float(delta_x_data)
    ppv_y = delta_y_plot/float(delta_y_data)

    coords = []
    coords_append = coords.append
    for x, y, err in data_points:

      if (y is not None) and (x is not None):
        x0 = (x - dr0)*ppv_x
        y0 = (y - dr1)*ppv_y
        x0 += pr0
        y0 += pr1
        
        if err is None:
          e1 = e2 = None
        else:
          e0 = err * ppv_y
          e1 = y0 + e0
          e2 = y0 - e0
        
        if y0 < pr3:
          continue
        if y0 > pr1:
          continue  
        
        #y0 = max(pr3, min(pr1, y0,) )
        
        coords_append((x0,y0,e1,e2))

    return coords
 
 
  def _graph_check_points(self, points):
    
    x_min = x_max = points[0][0]
    y_min = y_max = points[0][1]
    c = 0.0
    
    for i, point in enumerate(points):
      x, y = point[:2]
      
      x_min = min(x, x_min)
      y_min = min(y, y_min)
      x_max = max(x, x_max)
      y_max = max(y, y_max)
     
      
      if len(point) == 2:
        e = None
      else:
        e = point[2]

      if (not y) and (y != 0):
        y = None 

      if (not x) and (x != 0):
        x = None 

      else:
        if not isinstance(x, (float, int)):
          x = c
          c += 1.0

      points[i] = (x, y, e)   

    return x_min, x_max, y_min, y_max
  
  def _hex_to_rgb(self, hex_code):
    
    r = int(hex_code[1:3], 16)
    g = int(hex_code[3:5], 16)
    b = int(hex_code[5:7], 16)
    
    return r, g, b

  def _default_color_func(self, matrix, pos_color, neg_color, bg_color):
    
    n, m = matrix.shape
    color_matrix = np.zeros((n, m, 3), float)
    rgb0 = np.array(self._hex_to_rgb(bg_color), float)
    rgbP = np.array(self._hex_to_rgb(pos_color), float)
    rgbN = np.array(self._hex_to_rgb(neg_color), float)
    
    
    for i in range(n):
      for j in range(m):
        f = matrix[i,j]
        
        if f > 0:
          g = 1.0 - f
          color_matrix[i,j] = (g * rgb0) + (f * rgbP)
          
        elif f < 0:
          f *= -1
          g = 1.0 - f
          color_matrix[i,j] = (g * rgb0) + (f * rgbN)
        else:
          color_matrix[i,j] = rgb0

    color_matrix = np.clip(color_matrix, 0, 255)
    color_matrix = np.array(color_matrix, dtype=np.uint8)
  
    return color_matrix
    
  
  def density_matrix(self, matrix, width, height=None, x_grid=None, y_grid=None,
                     x_labels=None, y_labels=None, x_axis_label=None, y_axis_label=None, 
                     line_color='#000000', bg_color='#FFFFFF', grid_color='#808080',
                     pos_color='#0080FF', neg_color='#FF4000', color_func=None,
                     font=None, font_size=16, line_width=1, plot_offset=(50, 50),
                     value_range=None, scale_func=None):
    
    pad = font_size + font_size / 4
    
    if not height:
      height = width
    
    if not font:
      font = self.font  
    
    if not color_func:
      color_func = lambda x, p=pos_color, n=neg_color, b=bg_color: self._default_color_func(x, p, n, b)
    
    c_matrix = np.array(matrix)
    
    if value_range:
      a, b = value_range
      c_matrix = np.clip(c_matrix, a, b)
    
    if scale_func:
      c_matrix = scale_func(c_matrix)
    
    c_matrix /= max(c_matrix.max(), -c_matrix.min())
    
    n, m = c_matrix.shape
    
    x_box = width/float(m)
    y_box = height/float(n)
       
    x0, y0 = plot_offset
    
    x1, y1 = x0, y0
    
    if y_labels:
      x1 += pad
    
    if y_axis_label:
      x1 += pad
    
    x2 = x1 + width
    y2 = y1 + height
    
    c_matrix = color_func(c_matrix)
            
    self.rect((x1, y1, x2, y2), color=line_color, fill=bg_color)
    
    self.image(x1, y1, width, height, c_matrix)
            
    if x_grid:
      for val in x_grid:
        x = x1 + val * x_box
        self.line((x,y1,x,y2), color=grid_color, line_width=line_width)
        
    if y_grid:
      for val in y_grid:
        y = y1 + val * y_box
        self.line((x1,y,x2,y), color=grid_color, line_width=line_width)
    
    y3 = y2 + font_size
    x3 = x1 - font_size/2
    
    if x_labels:
      for i, val in enumerate(x_labels):
        if isinstance(val, (tuple, list)):
          x, t = val
        
        else:
          t = val  
          x = i
        
        x = min(float(m), max(0.0, x))
        
        x = x1 + +x_box/2.0 + x * x_box
    
        self.text(t, (x, y3), anchor='middle', size=font_size-2, bold=False, font=font, color=line_color, angle=None, vert_align=None)
      
      y3 += pad
       
    if y_labels:
      for i, val in enumerate(y_labels):
        if isinstance(val, (tuple, list)):
          y, t = val
        
        else:
          t = val  
          y = i
        
        y = min(float(n), max(0.0, y))
        y = y1 + y_box/2.0 + y * y_box
    
        self.text(t, (x3, y), anchor='middle', size=font_size-2, bold=False, font=font, color=line_color, angle=270, vert_align=None)
      
      x3 -= pad
    
    if x_axis_label:
      x = m/2.0
      x = x1 + x * x_box
      
      self.text(x_axis_label, (x, y3), anchor='middle', size=font_size, bold=False, font=font, color=line_color, angle=None, vert_align=None)
    
      y3 += pad 

    if y_axis_label:
      y = n/2.0
      y = y1 + y * y_box
    
      self.text(y_axis_label, (x3, y), anchor='middle', size=font_size, bold=False, font=font, color=line_color, angle=270, vert_align=None)
    
      x3 += pad 
         
         
         
  def graph(self, x, y, width, height, data_lists, x_label, y_label,
            names=None, colors=None,  graph_type=LINE_TYPE,
            symbols=None, line_widths=None, symbol_sizes=None,
            legend=False, title=None, x_labels=None, plot_offset=(100, 50),
            axis_color='black', bg_color='#F0F0F0', font=None, font_size=16, line_width=1,
            x_ticks=True, y_ticks=True, x_grid=False, y_grid=False):
    
    n_data = len(data_lists)
    
    if not names:
      names = ['Data %d' % (i+1) for i in range(n_data)]
    
    if not font:
      font = self.font
    
    if not colors:
      colors = DEFAULT_COLORS
    
    if not symbols:
      symbols = ['circle' for i in range(n_data)]
    
    if not line_widths:
      line_widths = [1 for i in range(n_data)]
    
    if not symbol_sizes:
      symbol_sizes = [2 for i in range(n_data)]
    
    # Data region, check santity
    x_min, x_max, y_min, y_max = self._graph_check_points(data_lists[0])
    
    for data in data_lists[1:]:
      x_min1, x_max1, y_min1, y_max1 = self._graph_check_points(data)
      x_min = min(x_min, x_min1)
      y_min = min(y_min, y_min1)
      x_max = max(x_max, x_max1)
      y_max = max(y_max, y_max1)
    
    data_region = (x_min, y_min, x_max, y_max)
    
    # Inner plot region
    
    off_x, off_y, = plot_offset 
    x0 = x + off_x
    y1 = y + off_y # Y increases down
    x1 = x0 + width - off_x
    y0 = y1 + height - off_y
  
    plot_region = (x0, y0, x1, y1)
    
    # Data region adjust
    
    if graph_type == BAR_TYPE:
      y_min = 0.0
      y_max *= 1.1
      x_min -= 0.5
      x_max += 0.5
    
    if x_min == x_max:
      x_max += 1.0
    
    if y_min == y_max:
      y_max += 1.0
    
    # Draw
    
    n_colors = len(colors)  

    self.rect(plot_region, axis_color, fill=bg_color)
    self.line((x0,y0,x1,y0), axis_color, line_width)
    self.line((x0,y0,x0,y1), axis_color, line_width)

    self._graph_draw_ticks(data_region, plot_region, x_label, y_label,
                           line_width, x_ticks, y_ticks,
                           x_grid, y_grid, x_labels)
    
    if legend:
      x2 = x0 + (1.05*(x1-x0))
      y2 = y0 + (0.90*(y1-y0))

      for i, data_set in enumerate(data_lists):
        name = names[i]
        color = colors[i % n_colors]
        symbol = symbols[i] if graph_type in (LINE_TYPE, SCATTER_TYPE) else 'square'
        self._graph_draw_symbols( [(x2, y2, None, None)], color, symbol, fill=color, symbol_size=7)
        self.text(name, (x2+font_size,y2), anchor='start', size=font_size)
        y2 += font_size
  
    for i, data_points in enumerate(data_lists):
      color = colors[i % n_colors]
      symbol = symbols[i]
      lw = line_widths[i]
      symbol_size = symbol_sizes[i]
      coords = self._graph_get_draw_coords(data_points, data_region, plot_region)
      
      if graph_type == LINE_TYPE:
        points = [point[:2] for point in coords]
        self.lines(points, color, lw)
        self._graph_draw_symbols(coords, color, symbol, symbol_size)
        
      elif graph_type == BAR_TYPE:
        self._graph_draw_boxes(coords, color, plot_region)
        self._graph_draw_symbols(coords, 'black', None, 8) # For error bars
             
      else: # SCATTER_TYPE
        self._graph_draw_symbols(coords, color, symbol, symbol_size, fill='#808080')

    if title:
      self.text(title, ((x0+x1)/2.0, y1/2.0))
 
  
  def _graph_draw_boxes(self, coords, color, plot_region):
  
    points = [point[:2] for point in coords]
    y1 = plot_region[1]
    x_list = [p[0] for p in points]
    x_list.sort()
    
    deltas = [x_list[i+1]-x_list[i] for i in range(len(x_list)-1)]
    deltas.sort()
    width = deltas[0]/3.0
    
    for x, y in points:
      x1 = x - width
      x2 = x + width
      y2 = y
      self.rect((x1,y1,x2,y2), color='#000000', fill=color)
  
  
  def _graph_draw_symbols(self, coords, color, symbol, symbol_size, fill='#808080'):
  
    radius = symbol_size/2.0
      
    if symbol == 'square':
      for x, y, yU, yL in coords:
        x0 = x-radius
        x1 = x+radius
        y0 = y-radius
        y1 = y+radius
    
        if yU is not None:
          self.line((x,y,x,yU), color=color)
          self.line((x0,yU,x1,yU), color=color)
 
        if yL is not None:
          self.line((x,y,x,yL), color=color)
          self.line((x0,yL,x1,yL), color=color)
       
        self.square((x,y), radius-1, color=color, fill=color)
   
    elif symbol is None:
      for x, y, yU, yL in coords:
        x0 = x-radius
        x1 = x+radius
        y0 = y-radius
        y1 = y+radius
    
        if yU is not None:
          self.line((x,y,x,yU), color=color)
          self.line((x0,yU,x1,yU), color=color)
 
        if yL is not None:
          self.line((x,y,x,yL), color=color)
          self.line((x0,yL,x1,yL), color=color)
         
    else:
      for x, y, yU, yL in coords:
        x0 = x-radius
        x1 = x+radius
        y0 = y-radius
        y1 = y+radius
        
        if yU is not None:
          self.line((x,y,x,yU), color=color)
          self.line((x0,yU,x1,yU), color=color)
 
        if yL is not None:
          self.line((x,y,x,yL), color=color)
          self.line((x0,yL,x1,yL), color=color)
         
        self.circle((x0, y0), radius, color=color, fill=color)
 
 
  def _graph_draw_ticks(self, data_region, plot_region, x_label, y_label,
                        line_width=1, x_ticks=True, y_ticks=True,
                        x_grid=False, y_grid=False, x_labels=None):
  
    def format_text(text):
 
      if isinstance(text, float):
        if text == 0:
          text = '0'
        elif abs(text) > 999999 or abs(text) < 0.01:
          text = '%5.2e' % text
        else:
          text = str(text)
          
      elif isinstance(text, int):
        text = str(text)
        
      if text and text[0:1] == '@':
        text = ''

      return text
    
    x0, y0, x1, y1 = plot_region
    delta_xplot = x1 - x0
    delta_yplot = y1 - y0
    delta_x_data = data_region[2] - data_region[0]
    delta_y_data = data_region[3] - data_region[1]
    size = 7
    
    if x_labels:
      sizes = [len(x) for x in x_labels]
      size = max(sizes)
  
    y_close = 40
    x_close = 140
    
    xs = x0 - 8
    ys = y0 + 8
 
    xt = xs - 4
    yt = ys + 16
    
    ppv_x = delta_xplot/float(delta_x_data)
    ppv_y = delta_yplot/float(delta_y_data)
  
    space_x_data = x_close/ppv_x
    space_y_data = y_close/ppv_y
    
    sci_x = '%e' % abs(space_x_data)
    sci_y = '%e' % abs(space_y_data)
    
    deci_x = int(sci_x[-3:])
    deci_y = int(sci_y[-3:])

    sig_d_x = int(sci_x[0])
    sig_d_y = int(sci_y[0])

    n_x = 10.0
    n_y = 10.0
    s_x = abs(sig_d_x-n_x)
    s_y = abs(sig_d_y-n_y)
    for n in (1.0, 2.0, 5.0):
      s = abs(sig_d_x-n)
      if s < s_x:
        s_x = s
        n_x = n

      s = abs(sig_d_y-n)
      if s < s_y:
        s_y = s
        n_y = n

    inc_x = (abs(space_x_data)/space_x_data) *  n_x * 10**(deci_x)
    inc_y = (abs(space_y_data)/space_y_data) * -n_y * 10**(deci_y)

    val_x = data_region[0] - (data_region[0] % inc_x )
    val_y = data_region[1] - (data_region[1] % inc_y )

    if x_ticks:
      for i in range(int(round(delta_x_data/inc_x))+2):
 
        tick_x = round(val_x,-deci_x)
        x = plot_region[0]+(tick_x - data_region[0])*ppv_x
        val_x += inc_x
          
        if x > plot_region[2]:
          continue
        if x < plot_region[0]:
          continue

        if x_labels and (i-1 < len(x_labels)):
          self.text(x_labels[i-1], (x, yt), anchor='middle')
        else:
          self.text(format_text(tick_x), (x, yt), anchor='middle')
 
        if x_grid:
          self.line((x, y0, x, y1), 'black', 1)
        
        self.line((x, y0, x, ys), 'black', 1)
    
    if y_ticks:
      for i in range(int(round(delta_y_data/inc_y))+2):
 
        tick_y = round(val_y,-deci_y)
        y = plot_region[1]+(tick_y - data_region[1])*ppv_y
        val_y += inc_y
 
        if y < plot_region[3]:
          continue
        if y > plot_region[1]:
          continue

        self.text(format_text(tick_y), (xt, y), anchor='end', vert_align='middle')
        
        if y_grid:
          self.line((x0, y, x1, y), 'black', 1)
        
        self.line((x0, y, xs, y), 'black', 1)

    if x_label:
      self.text(x_label, ((x0+x1)/2.0, yt+24))

    if y_label:
      self.text(y_label, (xt-65,((y0+y1)/2.0)), angle=-90)
   

  def pie_chart(self, x, y, height, values, texts=None, colors=None, line_color='black', small_val=0, line_width=1, box_size=16, pad=4):
    
    rad = height/2 - 4
    x0 = x + height/2
    y0 = y + height/2

    n = float(sum(values))
    nv = len(values)
    
    if not colors:
      colors = [colorsys.hsv_to_rgb(i/(nv+1.0), 0.7, 0.8) for i in range(nv)]
      colors = ['#%02X%02X%02X' % (int(r*255), int(g*255), int(b*255)) for r, g, b in colors]
    
    nc = len(colors)
    a0 = -pi/2
    other = 0.0
    c = 0
    
    for i, value in enumerate(values):
      if value/n < small_val:
        other += value
        continue  
        
      a1 = a0 + TAU * value/n
      self.segment(x0, y0, rad, a0, a1, line_color, colors[c%nc], line_width)
      a0 = a1
      c += 1
    
    if other:
      a1 = a0 + TAU * other/n
      self.segment(x0, y0, rad, a0, a1, line_color, colors[c%nc], line_width)
    
    x1 = x + height + 2*pad
    y1 = y + pad + box_size
    
    if texts:
    
      c = 0
      for i in range(nv):
        if values[i]/n < small_val:
          continue  
       
        self.rect((x1, y1, x1+box_size, y1-box_size), color=line_color, fill=colors[c%nc])
        self.text(texts[i], (x1+box_size+pad, y1-2), anchor='start', size=box_size)
        
        y1 += box_size + pad
        c += 1
    
      if other:
        self.rect((x1, y1, x1+box_size, y1-box_size), color=line_color, fill=colors[c%nc])
        self.text('Other', (x1+box_size+pad, y1-2), anchor='start', size=box_size)
    
  
  def table(self, x0, y0, width, data, header=True, text_anchors=None, col_formats=None,
            size=16, pad=2, font=None, main_color='black', line_color='#808080'):

    row_height = size + pad + pad
    if not font:
      font = self.font
   
    if isinstance(width, float) and width < 1.0:
      width *= 1000
    
    lines = []
    n_cols = len(data[0])
    
    if not col_formats:
      col_formats = ['%s' for x in range(n_cols)]

    if not text_anchors:
      text_anchors = ['start' for x in range(n_cols)]
    
    col_widths = defaultdict(int)
    for row in data:
       for col, t in enumerate(row):
         if t is None:
           continue
       
         t = col_formats[col] % t
         col_widths[col] = max(len(t), col_widths[col])
    
    n = float(sum(col_widths.values()))
    col_widths = [width * col_widths[i]/n for i in range(n_cols)]
    
    x, y, = x0, y0
    
    #self.group_start()   
    
    self.line((x, y, x+width, y), color=line_color, line_width=1)
    
    for i, row in enumerate(data):
      y += row_height
    
      if (i == 0) and header:
        self.line((x, y, x+width, y), color=line_color, line_width=1)
        bold = True
      else:
        bold = False  
      
      dx = 0
      for j, text in enumerate(row):
        if text is not None:
          text = col_formats[j] % text
          anchor = text_anchors[j]
          
          x1 = x + dx 
          
          if anchor == 'start':
            x1 += pad
            
          elif anchor == 'middle': 
          
            x1 += col_widths[j]/2
            
          else: # 'end'
            x1 += col_widths[j] - pad
                    
          self.text(text, (x1, y-pad), anchor, size, bold, font, color=main_color)
          
        dx += col_widths[j]
    
    y += pad
    y += pad
        
    self.line((x, y, x+width, y), color=line_color, line_width=1)
    
    #self.group_end()
    
    return width, y-y0

if __name__ == '__main__':


  svg_doc = SvgDocument()
  height = 1200
  
  """ 
  vals = [99, 44, 11]
  
  names = list('ABC')
  
  svg_doc.pie_chart(432, 804, 124, vals, names, colors=['#80C0FF','#FFFF00','#FF0000'], line_color='#808080')
  
  data_set = [(x, x**1.5) for x in range(100)]
  
  svg_doc.graph(432, 400, 500, 200, [data_set], 'Number', 'Value',
                names=None, colors=None,  graph_type='line',
                symbols=None, line_widths=None, symbol_sizes=None,
                legend=False, title=None, plot_offset=(100, 50))
  """
  
  xx = [x/10.0 for x in range(1,11)]
  yy = [x*x for x in xx]

  data = []
  for x in xx:
    data.append([])
    for y in yy:
      data[-1].append(x*y - 0.4)
  
  labels = list('ABCDEFGHIJ')
      
  svg_doc.density_matrix(data, 500, 500, x_grid=[5.0], y_grid=[5.0],
                   x_labels=labels, y_labels=labels, x_axis_label='X axis', y_axis_label='Y axis',
                   font=None, font_size=16, line_width=1, plot_offset=(50, 50),
                   value_range=None, scale_func=None)
   
  svg_doc.write_file('Test_chart.svg', height, height)

