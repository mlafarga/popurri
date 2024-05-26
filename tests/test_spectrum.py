from pathlib import Path

import numpy as np
import pandas as pd

from popurri import spectrum


# Functions to read spectra and extra data
# ----------------------------------------

def test_filprop_exists():
    assert Path(spectrum.filheader_kw_table).is_file(), "File does not exist"


def test_read_spectrum_properties():
    assert isinstance(spectrum.read_header_kw_table(), pd.DataFrame), "Output is not a DataFrame"


# Other utils
# -----------

# Example test
def test_dopplershift():
    result = spectrum.dopplershift(5000, 100000, rel=True)
    result_expected = 5001.668098731313
    assert np.isclose(result, result_expected, rtol=1e-05, atol=1e-08)


