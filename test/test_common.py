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


@pytest.mark.parametrize("x, n, expected", [
    ([1, 2, 3, 4, 5], 2, [(1, 2), (3, 4), (5, None)]),
    ([1, 2, 3, 4], 2, [(1, 2), (3, 4)]),
])
def test_nwise_longest(x, n, expected):
    assert list(nwise_longest(x, n))


@pytest.mark.parametrize("inputs, expected", [
    (['foo.r1.fq', 'foo.r2.fq'], 'foo.r1.r2.fq'),
    (['foo_bar.r1.fq', 'foo_baz.r2.fq'], 'foo_bar.baz.r1.r2.fq'),
    (['foo_bar.r1.fq', 'foo.baz.r2.fq'], 'foo_bar.baz.r1.r2.fq'),
    (['foo_bar_r1.fq', 'foo_baz.r2.fq'], 'foo_bar_baz.r1.r2.fq'),
    (['foo_bar.r1.fq', 'foo_baz.fq'], 'foo_bar.baz.r1.fq'),
    (['foo_bar.fq', 'foo_baz.r2.fq'], 'foo_bar.baz.r2.fq'),
])
def test_merge_file_names(inputs, expected):
    assert merge_file_names(inputs[0], inputs[1]) == expected


@pytest.mark.parametrize("inputs", [
    ('foo.r1.fq', 'foo.r2.fq.gz'),
    ('bar/foo.r1.fq', 'baz/foo.r2.fq'),
])
def test_merge_file_error(inputs):
    with pytest.raises(Exception):
        merge_file_names(inputs[0], inputs[1])
