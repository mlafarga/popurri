from pathlib import Path

import pandas as pd

from popurri import spectrograph


def test_filprop_exists():
    assert Path(spectrograph.filprop).is_file(), "File does not exist"

def test_read_spectrograph_properties():
    assert isinstance(spectrograph.read_spectrograph_properties(), pd.DataFrame), "Output is not a DataFrame"
