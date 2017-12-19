import gzip
import re
import os


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
  FILENAME_SPLIT_PATT = re.compile('[_\.]')

  # same dir, need non truncated name

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

  parts1 = FILENAME_SPLIT_PATT.split(file_root1)
  parts2 = FILENAME_SPLIT_PATT.split(file_root2)
  parts3 = []

  n1 = len(parts1)
  n2 = len(parts2)
  n = max(n1, n2)

  for i in range(n):

    if (i < n1) and (i < n2):
      a = parts1[i]
      b = parts2[i]

      parts3.append(a)
      if a != b:
        parts3.append(b)

    elif i < n1:
      parts3.append(parts1[i])
    else:
      parts3.append(parts2[i])

  file_root3 = sep.join(parts3)

  file_path3 = os.path.join(dir_name1, file_root3 + file_ext1)

  return file_path3
