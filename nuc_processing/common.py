import gzip
import re
import os
try:
  from itertools import zip_longest
except ImportError:
  # Python 2
  from itertools import izip_longest as zip_longest
from functools import partial


def open_file_r(file_path, buffering=-1):
  if file_path.endswith('.gz'):
    return gzip.open(file_path, 'rt')
  else:
    return open(file_path, 'rtU', buffering)


def strip_ext(name, ext):
  if name.endswith(ext):
    return name[:-len(ext)]
  return name


def nwise_longest(iterable, n=2, fillvalue=None):
  return zip_longest(*[iter(iterable)] * n, fillvalue=fillvalue)


def merge_file_names(file_path1, file_path2, sep='_'):
  SPLIT_PATT = re.compile('([_\.])')

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
  # Split on separators
  split_names = map(SPLIT_PATT.split, [file_root1, file_root2])
  # Pair segment and following separator
  split_names = map(partial(nwise_longest, n=2, fillvalue="."), split_names)

  for (a, sep_a), (b, sep_b) in zip_longest(*split_names, fillvalue=("", "")):
    if a is not None:
      parts.extend([a, sep_a])
    if a != b and b is not None:
      parts.extend([b, sep_b])

  return os.path.join(dir_name1, ''.join(parts) + file_ext1[1:])
