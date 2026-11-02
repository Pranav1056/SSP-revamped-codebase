import argparse
from ssp.astrometry import lspr
parser = argparse.ArgumentParser()
parser.add_argument("input_file")
parser.add_argument("x", type=float)
parser.add_argument("y", type=float)
args = parser.parse_args()
lspr(args.input_file, args.x, args.y)
