"""Original function names, retained for compatibility."""

from .orbital_elements import getMomentum, getAngle, getPerihelion, getA, getInclination, getEccentricity, getLongitudeOfAscendingNode, getTrueAnomaly, getArgOfPerihelion, getMeanAnomaly, getTimeOfLastPerihelion
from .kepler import newtonRaphson, fg, initialFg, scaleEqLagrange
from .coordinates import julian, convertToDecimalHours, convertToDecDegrees, getRaAndDec
from .astrometry import lspr
from .photometry import getCircle, aperturePhotometry
