import numpy as np
from ssp.kepler import newtonRaphson, initialFg

def test_kepler():
    e = newtonRaphson(0.7,0.2,1e-12)
    assert np.isclose(e-0.2*np.sin(e),0.7)

def test_initial_fg():
    assert initialFg(2,0)==(1,0)
