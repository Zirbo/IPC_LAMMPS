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
parser.add_argument('vEE', metavar='vEE', type=float, help='')
parser.add_argument('vEP', metavar='vEP', type=float, help='')
parser.add_argument('cPP', metavar='cPP', type=float, help='')
parser.add_argument('janus', metavar='janus', type=str, help='')
args = parser.parse_args()

ecc = args.ecc
delta = args.delta
vEE = args.vEE
vEP = args.vEP
cPP = args.cPP
janus = (args.janus == "janus");

HSradius = 0.5
HSdiameter = 1.0
bigRadius = HSradius + delta/2
patchRadius = bigRadius - ecc

def computeOmega(Ra, Rb, rab):
  """ overlap volumes: according to BKL paper, formula 18 """
  if ( rab > Ra+Rb ):
    return 0.
  elif ( rab <= fabs(Ra-Rb) ):
    return (min(Ra,Rb)/HSradius)**3;
  else:
    tempSum = (Ra**2-Rb**2)/(2.*rab)
    return (0.25/HSradius**3) * (
               (2.*Ra+tempSum+rab/2.)*(Ra-tempSum-rab/2.)**2 +
               (2.*Rb-tempSum+rab/2.)*(Rb+tempSum-rab/2.)**2 )


if __name__ == "__main__":

  # compute maximum overlap volumes
  fBB = computeOmega(bigRadius, bigRadius, HSdiameter)
  fBS = computeOmega(bigRadius, patchRadius, HSdiameter - ecc)
  fSS = computeOmega(patchRadius, patchRadius, HSdiameter - 2*ecc)
  
  # compute coefficients
  eBB = vEE
  eBS = vEP - eBB
  eSS = cPP - eBB - 2.*eBS
  cBB = eBB/fBB
  cBS = eBS/fBS
  cSS = eSS/fSS
  
  # create input file to generate the potentials
  outFile = open("inputfile_" + args.name + ".txt", 'w')
  if janus:
    outFile.write(str(cBB)[0:7].ljust(10) + str(cBS)[0:8].ljust(10) + str(0)[0:8].ljust(10) + "\n")
    outFile.write(str(cSS)[0:7].ljust(10) + str(0)[0:7].ljust(10) + str(0)[0:7].ljust(10) + "\n")
  else:
    outFile.write(str(cBB)[0:7].ljust(10) + str(cBS)[0:8].ljust(10) + str(cBS)[0:8].ljust(10) + "\n")
    outFile.write(str(cSS)[0:7].ljust(10) + str(cSS)[0:7].ljust(10) + str(cSS)[0:7].ljust(10) + "\n")
  outFile.write("1.0\n")
  outFile.write(str(ecc)[0:7].ljust(10) + str(patchRadius)[0:7].ljust(10) + "\n")
  outFile.write(str(ecc)[0:7].ljust(10) + str(patchRadius)[0:7].ljust(10) + "\n")
  outFile.write("500       15\n")
  outFile.write("1.0e-05   1.0e+06\n")

  outFile = open("MC_inputfile_" + args.name + ".txt", 'w')
  if janus:
    outFile.write(str(cBB)[0:7].ljust(10) + str(cBS)[0:8].ljust(10) + str(0)[0:8].ljust(10) + "\n")
    outFile.write(str(cSS)[0:7].ljust(10) + str(0)[0:7].ljust(10) + str(0)[0:7].ljust(10) + "\n")
  else:
    outFile.write(str(eBB)[0:7].ljust(10) + str(eBS)[0:8].ljust(10) + str(eBS)[0:8].ljust(10) + "\n")
    outFile.write(str(eSS)[0:7].ljust(10) + str(eSS)[0:7].ljust(10) + str(eSS)[0:7].ljust(10) + "\n")
  outFile.write("1.0\n")
  outFile.write(str(ecc)[0:7].ljust(10) + str(patchRadius)[0:7].ljust(10) + "\n")

