#! /usr/bin/python3

# writes a configuration with a plane structure of ipcs.

import argparse
from math import cos, sin, sqrt, pi, floor
from numpy.random import ranf

helpString = """Creates a LAMMPS starting configuration with a single wafer plane.\n
Suggested values for a cubic box: 14 12 1.2 1.2 12.4 12.4 0.22\n
Suggested values for an elongates box for gravity experiments: 14 12 1.2 sz 12.4 Lz 0.22\n"""

parser = argparse.ArgumentParser(description=helpString)
parser.add_argument('particlePerSideX', metavar='nPx', type=int, help='number of particles in the X side')
parser.add_argument('particlePerSideY', metavar='nPy', type=int, help='number of particles in the Y side')
parser.add_argument('zOrigin', metavar='z0', type=float, help='height of the plane')
parser.add_argument('density', metavar='rho', type=float, help='density')
parser.add_argument('boxSideX', metavar='Lx', type=float, help='size of the simulation box side base (x)')
parser.add_argument('boxSideY', metavar='Ly', type=float, help='size of the simulation box side base (y)')
parser.add_argument('boxSideZ', metavar='Lz', type=float, help='height of the simulation box side (z)')
parser.add_argument('ecc', metavar='e', type=float, help='eccentricity of the IPCs')
args = parser.parse_args()
print(args)

outputFile = open('startingstate_manner_plus_FCC.txt','w')

Lx = args.boxSideX
Ly = args.boxSideY
Lz = args.boxSideZ
zOrigin = args.zOrigin
if zOrigin < 0.5:
    zOrigin = 0.5
density = args.density
ecc = args.ecc
nWaferX = args.particlePerSideX
nWaferY = args.particlePerSideY
nIPCs = int(density*Lx*Ly*Lz)
totalWaferParticles = nWaferX*nWaferY
totalFluidParticles = nIPCs - totalWaferParticles
spacing = 1.6
LR = spacing*.5
nFluidX = int(Lx/spacing)
nFluidY = int(Ly/spacing)
nFluidZ = int((Lz-zOrigin)/spacing)


def absolutePBCx(x):
    return x - Lx*floor(x/Lx)

def absolutePBCy(y):
    return y - Ly*floor(y/Ly)

def absolutePBCz(z):
    return z - Lz*floor(z/Lz)

alpha = .45*pi
beta = .93*pi
gamma = 0.25*pi
cos30 = sqrt(3)*.5
p = [ [ ecc*cos(alpha), ecc*sin(alpha), 0. ] ,
      [ ecc*cos(beta),  ecc*sin(beta),  0. ] ,
      [ ecc*cos(gamma), ecc*sin(gamma), 0. ] ]


outputFile.write("# 3D starting configuration for LAMMPS created with a script available at\n")
outputFile.write("# https://github.com/Zirbo/OSPC_LAMMPS\n")
outputFile.write("# The plane particles are from 1 to " + str(nWaferX*nWaferY))

outputFile.write("\n")
outputFile.write("\n" + str(3*nIPCs).rjust(16) + " atoms")
outputFile.write("\n" + str(2*nIPCs).rjust(16) + " bonds")
outputFile.write("\n" + str(  nIPCs).rjust(16) + " angles")
outputFile.write("\n")

outputFile.write("\n" + str(2).rjust(16) + " atom types")
outputFile.write("\n" + str(1).rjust(16) + " bond types")
outputFile.write("\n" + str(1).rjust(16) + " angle types")
outputFile.write("\n")

outputFile.write("\n" + '{:3.8f}'.format(0.0).rjust(16) +
                        '{:3.8f}'.format(Lx).rjust(16) + "     xlo xhi")
outputFile.write("\n" + '{:3.8f}'.format(0.0).rjust(16) +
                        '{:3.8f}'.format(Ly).rjust(16) + "     ylo yhi")
outputFile.write("\n" + '{:3.8f}'.format(0.0).rjust(16) +
                        '{:3.8f}'.format(Lz).rjust(16) + "     zlo zhi")

outputFile.write("\n")
outputFile.write("\nMasses")
outputFile.write("\n#  atomtype, mass")
outputFile.write("\n" + str(1).rjust(10) + str(2.0).rjust(10))
outputFile.write("\n" + str(2).rjust(10) + str(0.5).rjust(10))

outputFile.write("\n")
outputFile.write("\nAtoms")
outputFile.write("\n#   atom-ID    mol-ID   atom-type    charge    x               y               z")

# wafer layer
waferParticles = 0
z = zOrigin
for ix in range(nWaferX):
    x = (0.6 + 1.0000000000001*cos30*ix)
    for iy in range(nWaferY):
        waferParticles += 1
        atomNumber = (waferParticles - 1)*3 + 1
        # ipc center
        j = 0 if (iy + (int((ix + 1)/2))%2)%2==0 else 1
        y = 0.6 + ( (.5 + 1.0000000000001*iy) if ix%2==0 else (1.0000000000001*iy) )
        outputFile.write("\n" + str(atomNumber).rjust(10) +
              str(waferParticles).rjust(10) +
              str(1).rjust(10) +
              str(-1.).rjust(10) +
             '{:3.8f}'.format(x).rjust(16) +
             '{:3.8f}'.format(y).rjust(16) +
             '{:3.8f}'.format(z).rjust(16) )
        # first patch
        px = x + p[j][0];    px = absolutePBCx(px)
        py = y + p[j][1];    py = absolutePBCy(py)
        pz = z + p[j][2];    pz = absolutePBCz(pz)
        atomNumber += 1
        outputFile.write("\n" + str(atomNumber).rjust(10) +
              str(waferParticles).rjust(10) +
              str(2).rjust(10) +
              str(0.5).rjust(10) +
             '{:3.8f}'.format(px).rjust(16) +
             '{:3.8f}'.format(py).rjust(16) +
             '{:3.8f}'.format(pz).rjust(16) )
        # second patch
        px = x - p[j][0];    px = absolutePBCx(px)
        py = y - p[j][1];    py = absolutePBCy(py)
        pz = z - p[j][2];    pz = absolutePBCz(pz)
        atomNumber += 1
        outputFile.write("\n" + str(atomNumber).rjust(10) +
              str(waferParticles).rjust(10) +
              str(2).rjust(10) +
              str(0.5).rjust(10) +
             '{:3.8f}'.format(px).rjust(16) +
             '{:3.8f}'.format(py).rjust(16) +
             '{:3.8f}'.format(pz).rjust(16) )

# square lattice above the plane
lattice_roots=( (0., 0., 0., 0),
                (LR, LR, 0., 1),
                (LR, 0., LR, 2),
                (0., LR, LR, 3) )

fluidParticles = 0
px = ecc/(3**.5)
py = ecc/(3**.5)
pz = ecc/(3**.5)
zStart = zOrigin + 1.2
for iz in range(nFluidZ):
    for iy in range(nFluidY):
        for ix in range(nFluidX):
            for r in lattice_roots:
                # ipc center
                fluidParticles += 1
                atomNumber = waferParticles*3 + (fluidParticles - 1)*3 + 1
                x = 0.5    + r[0] + ix*spacing + 0.003*ranf()
                y = 0.5    + r[1] + iy*spacing + 0.003*ranf()
                z = zStart + r[2] + iz*spacing + 0.003*ranf()
                outputFile.write("\n" + str(atomNumber).rjust(10) +
                      str(waferParticles + fluidParticles).rjust(10) +
                      str(1).rjust(10) +
                      str(-1.).rjust(10) +
                     '{:3.8f}'.format(x).rjust(16) +
                     '{:3.8f}'.format(y).rjust(16) +
                     '{:3.8f}'.format(z).rjust(16) )
                # create random unit vector
                px = ranf() 
                py = ranf()
                pz = ranf()
                mod = ecc/ ( (px**2 + py**2 + pz**2)**(0.5) )
                px *= mod
                py *= mod
                pz *= mod
                # patches
                outputFile.write("\n" + str(atomNumber + 1).rjust(10) +
                      str(waferParticles + fluidParticles).rjust(10) +
                      str(2).rjust(10) +
                      str(0.5).rjust(10) +
                     '{:3.8f}'.format(x + px).rjust(16) +
                     '{:3.8f}'.format(y + py).rjust(16) +
                     '{:3.8f}'.format(z + pz).rjust(16) )
                outputFile.write("\n" + str(atomNumber + 2).rjust(10) +
                      str(waferParticles + fluidParticles).rjust(10) +
                      str(2).rjust(10) +
                      str(0.5).rjust(10) +
                     '{:3.8f}'.format(x - px).rjust(16) +
                     '{:3.8f}'.format(y - py).rjust(16) +
                     '{:3.8f}'.format(z - pz).rjust(16) )
                if ( fluidParticles == totalFluidParticles):
                   break
            if ( fluidParticles == totalFluidParticles):
               break
        if ( fluidParticles == totalFluidParticles):
           break
    if ( fluidParticles == totalFluidParticles):
       break


outputFile.write("\n")
outputFile.write("\nBonds")
outputFile.write("\n#  ID bond-type atom-1 atom-2")
for i in range(nIPCs):
    IDcenter = 3*i + 1
    IDpatch1 = 3*i + 2
    IDpatch2 = 3*i + 3
    outputFile.write("\n" + str(2*i+1).rjust(10) + str(1).rjust(10) +
                            str(IDcenter).rjust(10) + str(IDpatch1).rjust(10) )
    outputFile.write("\n" + str(2*i+2).rjust(10) + str(1).rjust(10) +
                            str(IDcenter).rjust(10) + str(IDpatch2).rjust(10) )

outputFile.write("\n")
outputFile.write("\nAngles")
outputFile.write("\n#  ID    angle-type atom-1 atom-2 atom-3  (atom-2 is the center atom in angle)")
for i in range(nIPCs):
    IDcenter = 3*i + 1
    IDpatch1 = 3*i + 2
    IDpatch2 = 3*i + 3
    outputFile.write("\n" + str(i+1).rjust(10) + str(1).rjust(10) +
                            str(IDpatch1).rjust(10) + str(IDcenter).rjust(10) +
                            str(IDpatch2).rjust(10) )
outputFile.write("\n")
