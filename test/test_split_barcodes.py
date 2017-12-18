from nuc_processing.splitFastqBarcodes import *


def test_split(tmpdir):
    main(
        ['test/fixtures/sample.s_1.r_1.fq', 'test/fixtures/sample.s_1.r_2.fq', str(tmpdir)]
    )
    assert len(tmpdir.listdir()) == 12
    for outfile in tmpdir.listdir():
        lines = outfile.readlines()
        assert len(lines) >= 4
        assert len(lines) % 4 == 0
