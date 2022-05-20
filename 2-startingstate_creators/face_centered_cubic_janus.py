#! /usr/bin/python3

# writes a configuration with a face-centered-cubic lattice structure of janus OSPCs

import argparse
from math import cos, sin, sqrt, pi
from numpy.random import ranf

parser = argparse.ArgumentParser(description='Creates a LAMMPS starting configuration.')
parser.add_argument('nOSPCsSide', metavar='N', type=int, help='number of OSPCs: 4*N**2')
parser.add_argument('density', metavar='d', type=float, help='density')
parser.add_argument('ecc', metavar='e', type=float, help='eccentricity off the OSPCs')
args = parser.parse_args()

args.nOSPCsPlane = args.nOSPCsSide**2
args.nOSPCsCube  = args.nOSPCsSide**3
args.nOSPCs = 4*args.nOSPCsCube
args.side = (args.nOSPCs/args.density)**(1./3.)
args.sideStep = args.side/args.nOSPCsSide

print(args)

filename = "startingstate_FCC_1p_{}_{}_{}.txt".format(args.nOSPCs, args.density, args.ecc)
outputFile = open(filename,'w')

lattice_roots=( (0.0, 0.0, 0.0, 0),
                (0.5, 0.5, 0.0, 1),
                (0.5, 0.0, 0.5, 2),
                (0.0, 0.5, 0.5, 3) )



outputFile.write("# starting configuration for LAMMPS generated with a script available at\n")
outputFile.write("# https://github.com/Zirbo/OSPC_LAMMPS")
outputFile.write("\n# FCC crystal 1-patch DENSITY: {}, ecc: {}\n".format(args.density, args.ecc))
outputFile.write("\n" + str(2*args.nOSPCs).rjust(16) + " atoms")
outputFile.write("\n" + str(1*args.nOSPCs).rjust(16) + " bonds")
outputFile.write("\n")

outputFile.write("\n" + str(2).rjust(16) + " atom types")
outputFile.write("\n" + str(1).rjust(16) + " bond types")
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
outputFile.write("\n" + str(1).rjust(10) + str(2.9).rjust(10))
outputFile.write("\n" + str(2).rjust(10) + str(0.1).rjust(10))

outputFile.write("\n")
outputFile.write("\nAtoms")
outputFile.write("\n#  atom-ID mol-ID atom-type charge    x               y               z")
for r in lattice_roots:
    for iz in range(0, args.nOSPCsSide):
        for iy in range(0, args.nOSPCsSide):
            for ix in range(0, args.nOSPCsSide):
                # ospc center
                i = 1 + 2*(ix + iy*args.nOSPCsSide + iz*args.nOSPCsPlane + args.nOSPCsCube*r[3])
                ospcNum = int(i/2) + 1
                x = (0.5 + r[0] + ix + 0.03*ranf())*args.sideStep
                y = (0.5 + r[1] + iy + 0.03*ranf())*args.sideStep
                z = (0.5 + r[2] + iz + 0.03*ranf())*args.sideStep
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
                      str(1.0).rjust(10) +
                     '{:3.8f}'.format(x + args.ecc*px).rjust(16) +
                     '{:3.8f}'.format(y + args.ecc*py).rjust(16) +
                     '{:3.8f}'.format(z + args.ecc*pz).rjust(16) )

outputFile.write("\n")
outputFile.write("\nBonds")
outputFile.write("\n#  ID bond-type atom-1 atom-2")
for i in range(args.nOSPCs):
    IDcenter = 2*i + 1
    IDpatch = 2*i + 2
    outputFile.write("\n" + str(i+1).rjust(10) + str(1).rjust(10) +
                            str(IDcenter).rjust(10) + str(IDpatch).rjust(10) )

outputFile.write("\n")
