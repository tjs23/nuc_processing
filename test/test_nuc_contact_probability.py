from nuc_processing.NucContactMap import *

import pytest


def test_help():
    with pytest.raises(SystemExit) as e:
        main(["--help"])
        assert e.value.code == 0
