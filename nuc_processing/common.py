import gzip
import re
import os
try:
  from itertools import zip_longest
except ImportError:
  # Python 2
  from itertools import izip_longest as zip_longest


def open_file_r(file_path, buffering=-1):
  if file_path.endswith('.gz'):
    return gzip.open(file_path, 'rt')
  else:
    return open(file_path, 'rtU', buffering)


def strip_ext(name, ext):
  if name.endswith(ext):
    return name[:-len(ext)]
  return name


def merge_file_names(file_path1, file_path2, sep='_'):
  SPLIT_PATT = re.compile('[_\.]')

  dir_name1, file_name1 = os.path.split(file_path1)
  dir_name2, file_name2 = os.path.split(file_path2)

  if dir_name1 != dir_name2:
    msg = 'Attempt to merge file names for file from different directories'
    raise Exception(msg)

  file_root1, file_ext1 = os.path.splitext(file_name1)
  file_root2, file_ext2 = os.path.splitext(file_name2)

  if file_ext1 != file_ext2:
    msg = 'Attempt to merge file names with different file extensions'
    raise Exception(msg)

  parts = []
  for a, b in zip_longest(*map(SPLIT_PATT.split, [file_root1, file_root2])):
    parts.append(a)
    if a != b:
      parts.append(b)

  return os.path.join(dir_name1, sep.join(parts) + file_ext1)
