# OSPC for LAMMPS

Implementation of Overlapping Sphere Patchy Colloids for LAMMPS.

We provide all that is necessary to simulate OSPCs with LAMMPS:
sample .in files and startingstates, a program to print and tabulate the
OSPC potential, and generators of startingstates.
Anybody can start simulating IPCs in LAMMPS in few minutes.
For any help, feel free to contact us here on GitHub.

We have the following directories:

## 1-make_potentials
Contains programs that compute the potentials in LAMMPS format, and some helpers.

### from_contact_values_to_potentials.sh
(requires a c++14-compatible version of g++ and uses Eigen library,
which will be downloaded automatically)

Requires as input (modify the file to supply the inputs!)
the geometric parameters and the desired contact values of the potential
in the three EE, EP1 and P1P1 orientations
(plus EP2, P1P2, P2P2 for asymmetric models)
or BB, BP and PP for Janus,
and prints out the potential in an output directory,
named `target_<modelname>_<symmetry>_contact`.

The directory contains a lammpspot dir that contains your potential.
In addition there will be three directories:
- lammpspot_angular_plots
- lammpspot_radial_plots
- potential_on_path
whose content can be plotted using gnuplot or xmgrace, so you can
quickly visualize your potential. See the reference paper [3] for an
explanation of the orientations.

Finally there is also a file, inputfile.dat,
that recaps your input values.

### from_epsilon_to_potentials.sh
Requires as input (modify the file to supply the inputs!)
the geometric parameters and the epsilons from the mapping (see the
reference papers!) and prints out the potential in the output directory
named `target_<modelname>_<symmetry>_epsilons`.

It is a wrapper to the same program used by from_contact_values_to_potentials.sh.
Outputs and requirements are exactly the same.


### compute_ipc_geometry.py
(requires python3)

In case you want to follow the IPC geometry, this python script helps
you determine eccentricity and patch radius from patch amplitude and
interaction range, or viceversa. Check the inline help with -h.

## 2-startingstate_creators
Contains Python3 scripts that can be used to generate FCC startingstates
for Janus and two-patch OSPCs. There's also one cubic lattice for two-patch only.
All of them have an inline help, so run them with the -h flag to get explanations.

## 3-lammps_inputfiles
Each subdirectory of it contains a sample .in file and a sample
startingstate.
To run run a simulation you have to:

- copy th .in file to the directory where you want to run
- copy in that directory the potential and startingstate that you want to use
- adapt the .in file with paths to startingstate, potentials and eccentricities.
- run LAMMPS:
```
$ mpirun -np <num_cpus> </lammps/exe/path> -in run_*.in -var seed $RANDOM
```
where the seed is used to initialize the velocities (`$RANDOM` is a bash command).
Your LAMMPS distribution needs to be compiled with the MOLECULE package.

There is also a sample tk script that you can load in VMD (Visual Molecular
Dynamics) to visualize your trajectory:
```
$ vmd -e vmdscript.tk </path/to/trajectory.lammpstrj>
```

## 4-postprocess
!!! currently WIP !!!

 - does not support Janus
 - does not support Asymmetric OSPCs
 - works with any symmetric OSPC, although the code still refers to IPCs

Run the build.sh script to obtain the binary;
a c++11-compatible version of g++ and CMake 3.5 are required.

Run the obtained binary without arguments to see what arguments are required.
You need a valid trajectory, the directory with the potentials used to run the
simulation, and an inputfile with eccentricity and patch radius.

Will be extended to all OSPCs as soon as possible :)

## 9-advanced
Contains other startingstate generators and visualization tools
that we have been using in our research.
You are welcome to peep, but they are not all well documented,
so use them at your own risk ;)


## Have fun!
We hope that you will have as much fun playing with OSPCs as we did :D
S.F. & E.B.

## References
[1] "Inverse patchy colloids: from microscopic description to mesoscopic coarse-graining"
E. Bianchi, G. Kahl, and C. N. Likos - Soft Matter 7, 8313 (2011) 

[2] "Molecular dynamics simulations of inverse patchy colloids"
S. Ferrari, G. Kahl, E. Bianchi - Eur. Phys. J. E 41, 43 (2018) 

[3] ""
S. Ferrari, E. Locatelli, E. Bianchi - in prep.
