from nuc_processing.common import *

import pytest


def test_open_file_r(tmpdir):
    contents = "foo\nbar"
    uncompressed = str(tmpdir.join("test"))
    compressed = str(tmpdir.join("test.gz"))

    with open(uncompressed, 'w') as f:
        f.write(contents)
    with gzip.open(compressed, 'w') as f:
        f.write(contents.encode())

    for path in [uncompressed, compressed]:
        with open_file_r(path) as f:
            assert f.read() == contents


@pytest.mark.parametrize("x, ext, expected", [
    ('foo.gz', '.gz', 'foo'),
    ('foo.gz', '.gzip', 'foo.gz'),
])
def test_strip_ext(x, ext, expected):
    assert strip_ext(x, ext) == expected
