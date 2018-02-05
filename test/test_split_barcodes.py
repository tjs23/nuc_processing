from nuc_processing.splitFastqBarcodes import *
import gzip
import sys

if sys.version_info.major < 3:
    from pathlib2 import Path
else:
    from pathlib import Path

fixtures = list(map(Path, ['test/fixtures/sample.s_1.r_1.fq', 'test/fixtures/sample.s_1.r_2.fq']))


def test_split(tmpdir):
    infiles = fixtures
    main(list(map(str, infiles)) + [str(tmpdir)])
    assert len(tmpdir.listdir()) == 12
    for outfile in tmpdir.listdir():
        lines = outfile.readlines()
        assert len(lines) >= 4
        assert len(lines) % 4 == 0


def test_split_compressed(tmpdir):
    infiles = [Path(str(tmpdir)) / o.with_suffix('.fq.gz').name for o in fixtures]
    for original, infile in zip(fixtures, infiles):
        with original.open("rb") as original_f, gzip.open(str(infile), "wb") as in_f:
            in_f.write(original_f.read())
    main(list(map(str, infiles)) + [str(tmpdir)])
    for infile in infiles:
        infile.unlink()

    assert len(tmpdir.listdir()) == 12
    for outfile in tmpdir.listdir():
        lines = outfile.readlines()
        assert len(lines) >= 4
        assert len(lines) % 4 == 0
