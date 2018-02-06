from nuc_processing.NucProcess import *
import pytest
import sys

if sys.version_info.major < 3:
    from pathlib2 import Path
else:
    from pathlib import Path

reference = list(map(Path, ['test/fixtures/chr18.fa', 'test/fixtures/chr19.fa']))
reads = list(map(Path, ['test/fixtures/try_1.r_1_AAA.fq', 'test/fixtures/try_1.r_2_AAA.fq']))


def test_help():
    with pytest.raises(SystemExit) as e:
        main(["--help"])
        assert e.value.code == 0


def test_process(tmpdir):
    tmpdir = Path(str(tmpdir))
    main(
        ["-i"] + list(map(str, reads)) + ['-f'] + list(map(str, reference)) +
        ['-g', str(tmpdir / "ref")] + ['-o', str(tmpdir)]
    )
