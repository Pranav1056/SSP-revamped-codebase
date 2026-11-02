import numpy as np
from ssp.coordinates import julian, convertToDecDegrees

def test_julian():
    assert np.isclose(julian(False,2000,1,1,"12:00:00"),2451545)

def test_ra():
    assert convertToDecDegrees("12:00:00",False,True)==180
