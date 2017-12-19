import gzip


def open_file_r(file_path, buffering=-1):
  if file_path.endswith('.gz'):
    return gzip.open(file_path, 'rt')
  else:
    return open(file_path, 'rtU', buffering)


def strip_ext(name, ext):
  if name.endswith(ext):
    return name[:-len(ext)]
  return name
