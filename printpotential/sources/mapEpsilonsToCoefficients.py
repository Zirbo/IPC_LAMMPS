#! /usr/bin/python3


import argparse
from math import sqrt, pi, fabs, cos
from math import sin as sen

helpString = """Computes e_ij coefficients from the desired contact values.

Sample calls:
45n  0.2 0.22 0.1 -1 4       (45n)
"""

parser = argparse.ArgumentParser(description=helpString)
parser.add_argument('name', metavar='m', type=str, help='model name')
parser.add_argument('delta', metavar='d', type=float, help='interaction range minus diameter')
parser.add_argument('ecc', metavar='e', type=float, help='eccentricity')
parser.add_argument('eEE', metavar='eEE', type=float, help='')
parser.add_argument('eEP', metavar='eEP', type=float, help='')
parser.add_argument('ePP', metavar='ePP', type=float, help='')
parser.add_argument('emin', metavar='emin', type=float, help='')
args = parser.parse_args()

ecc = args.ecc
delta = args.delta
eEE = args.eEE
eEP = args.eEP
ePP = args.ePP
emin = args.emin

HSradius = 0.5
HSdiameter = 1.0
bigRadius = HSradius + delta/2
patchRadius = bigRadius - ecc


if __name__ == "__main__":

  # create input file to generate the potentials
  outFile = open("inputfile_" + args.name + ".txt", 'w')
  outFile.write(str(eEE)[0:7].ljust(10) + str(eEP)[0:8].ljust(10) + str(eEP)[0:8].ljust(10) + "\n")
  outFile.write(str(ePP)[0:7].ljust(10) + str(ePP)[0:7].ljust(10) + str(ePP)[0:7].ljust(10) + "\n")
  outFile.write(str(emin)[0:7].ljust(10) + "\n")
  outFile.write(str(ecc)[0:7].ljust(10) + str(patchRadius)[0:7].ljust(10) + "\n")
  outFile.write(str(ecc)[0:7].ljust(10) + str(patchRadius)[0:7].ljust(10) + "\n")
  outFile.write("500       15\n")
  outFile.write("1.0e-05   1.0e+06\n")

