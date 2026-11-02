import numpy as np
from ssp.orbital_elements import getMomentum, getA, getEccentricity, getAngle

def test_momentum():
    np.testing.assert_array_equal(getMomentum([1,0,0],[0,1,0]),[0,0,1])

def test_circular_orbit():
    assert np.isclose(getA([1,0,0],[0,1,0]),1)
    assert np.isclose(getEccentricity([0,0,1],1),0)

def test_angle():
    assert getAngle(0,1)==90
