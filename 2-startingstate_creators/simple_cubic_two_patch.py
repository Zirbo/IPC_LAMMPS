#! /usr/bin/python3

# writes a configuration with a cubic lattice structure of OSPCs

import argparse
from math import cos, sin, sqrt, pi
from numpy.random import ranf

parser = argparse.ArgumentParser(description='Creates a LAMMPS starting configuration.')
parser.add_argument('nOSPCsSide', metavar='N', type=int, help='cubic root of the number of particles')
parser.add_argument('density', metavar='d', type=float, help='density')
parser.add_argument('ecc1', metavar='e1', type=float, help='eccentricity of the first patch')
parser.add_argument('ecc2', metavar='e2', type=float, help='eccentricity of the second patch')
args = parser.parse_args()

outputFile = open('startingstate.txt','w')
args.nOSPCsPlane = args.nOSPCsSide**2
args.nOSPCs = args.nOSPCsSide**3
args.side = (args.nOSPCs/args.density)**(1./3.)
args.sideStep = args.side/args.nOSPCsSide

print(args)


outputFile.write("# starting configuration for LAMMPS generated with a script available at\n")
outputFile.write("# https://github.com/Zirbo/OSPC_LAMMPS")
outputFile.write("\n")
outputFile.write("\n" + str(3*args.nOSPCs).rjust(16) + " atoms")
outputFile.write("\n" + str(2*args.nOSPCs).rjust(16) + " bonds")
outputFile.write("\n" + str(  args.nOSPCs).rjust(16) + " angles")
outputFile.write("\n")

outputFile.write("\n" + str(3).rjust(16) + " atom types")
outputFile.write("\n" + str(2).rjust(16) + " bond types")
outputFile.write("\n" + str(1).rjust(16) + " angle types")
outputFile.write("\n")

outputFile.write("\n" + '{:3.8f}'.format(0.0).rjust(16) +
                        '{:3.8f}'.format(args.side).rjust(16) + "     xlo xhi")
outputFile.write("\n" + '{:3.8f}'.format(0.0).rjust(16) +
                        '{:3.8f}'.format(args.side).rjust(16) + "     ylo yhi")
outputFile.write("\n" + '{:3.8f}'.format(0.0).rjust(16) +
                        '{:3.8f}'.format(args.side).rjust(16) + "     zlo zhi")

outputFile.write("\n")
outputFile.write("\nMasses")
outputFile.write("\n#  atomtype, mass")
outputFile.write("\n" + str(1).rjust(10) + str(2.0).rjust(10))
outputFile.write("\n" + str(2).rjust(10) + str(0.5).rjust(10))
outputFile.write("\n" + str(3).rjust(10) + str(0.5).rjust(10))

outputFile.write("\n")
outputFile.write("\nAtoms")
outputFile.write("\n#  atom-ID mol-ID atom-type charge    x               y               z")
for iz in range(0, args.nOSPCsSide):
    for iy in range(0, args.nOSPCsSide):
        for ix in range(0, args.nOSPCsSide):
            # ospc center
            i = 1 + 3*(ix + iy*args.nOSPCsSide + iz*args.nOSPCsPlane)
            ospcNum = int(i/3) + 1
            x = (0.5 + ix + 0.1*ranf())*args.sideStep
            y = (0.5 + iy + 0.1*ranf())*args.sideStep
            z = (0.5 + iz + 0.1*ranf())*args.sideStep
            outputFile.write("\n" + str(i).rjust(10) +
                  str(ospcNum).rjust(10) +
                  str(1).rjust(10) +
                  str(-1.).rjust(10) +
                 '{:3.8f}'.format(x).rjust(16) +
                 '{:3.8f}'.format(y).rjust(16) +
                 '{:3.8f}'.format(z).rjust(16) )
            # create random unit vector
            px = ranf() 
            py = ranf()
            pz = ranf()
            mod = (px**2 + py**2 + pz**2)**(0.5)
            px /= mod
            py /= mod
            pz /= mod
            # patches
            outputFile.write("\n" + str(i+1).rjust(10) +
                  str(ospcNum).rjust(10) +
                  str(2).rjust(10) +
                  str(0.5).rjust(10) +
                 '{:3.8f}'.format(x + args.ecc1*px).rjust(16) +
                 '{:3.8f}'.format(y + args.ecc1*py).rjust(16) +
                 '{:3.8f}'.format(z + args.ecc1*pz).rjust(16) )
            outputFile.write("\n" + str(i+2).rjust(10) +
                  str(ospcNum).rjust(10) +
                  str(3).rjust(10) +
                  str(0.5).rjust(10) +
                 '{:3.8f}'.format(x - args.ecc2*px).rjust(16) +
                 '{:3.8f}'.format(y - args.ecc2*py).rjust(16) +
                 '{:3.8f}'.format(z - args.ecc2*pz).rjust(16) )

outputFile.write("\n")
outputFile.write("\nBonds")
outputFile.write("\n#  ID bond-type atom-1 atom-2")
for i in range(args.nOSPCs):
    IDcenter = 3*i + 1
    IDpatch1 = 3*i + 2
    IDpatch2 = 3*i + 3
    outputFile.write("\n" + str(2*i+1).rjust(10) + str(1).rjust(10) +
                            str(IDcenter).rjust(10) + str(IDpatch1).rjust(10) )
    outputFile.write("\n" + str(2*i+2).rjust(10) + str(2).rjust(10) +
                            str(IDcenter).rjust(10) + str(IDpatch2).rjust(10) )

outputFile.write("\n")
outputFile.write("\nAngles")
outputFile.write("\n#  ID    angle-type atom-1 atom-2 atom-3  (atom-2 is the center atom in angle)")
for i in range(args.nOSPCs):
    IDcenter = 3*i + 1
    IDpatch1 = 3*i + 2
    IDpatch2 = 3*i + 3
    outputFile.write("\n" + str(i+1).rjust(10) + str(1).rjust(10) +
                            str(IDpatch1).rjust(10) + str(IDcenter).rjust(10) +
                            str(IDpatch2).rjust(10) )
outputFile.write("\n")
