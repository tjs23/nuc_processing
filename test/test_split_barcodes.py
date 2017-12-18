from nuc_processing.splitFastqBarcodes import *

import pytest


@pytest.mark.parametrize("infiles", [
    ['test/fixtures/sample.s_1.r_1.fq', 'test/fixtures/sample.s_1.r_2.fq'],
    ['test/fixtures/sample.s_1.r_1.fq.gz', 'test/fixtures/sample.s_1.r_2.fq.gz'],
])
def test_split(tmpdir, infiles):
    main(infiles + [str(tmpdir)])
    assert len(tmpdir.listdir()) == 12
    for outfile in tmpdir.listdir():
        lines = outfile.readlines()
        assert len(lines) >= 4
        assert len(lines) % 4 == 0
