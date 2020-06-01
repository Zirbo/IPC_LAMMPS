IPC for LAMMPS
Implementation of Inverse Patchy Colloids for LAMMPS.

We provide sample .in files and startingstates for LAMMPS, and
a program to print a tabulate of the IPC potential, so that
anybody interested can start simulating IPCs in LAMMPS in few minutes :)

For any help, feel free to contact us!

LAST UPDATED: 1 June 2020

We have 4 directories.

- inputfiles: each subdirectory of it contains a sample .infile and the
  startingstates that it requires.
   1 - Copy them to the directory where you want to run.
   2 - Copy in that directory the potential that you want to use,
         (see printpotential below!)
     and fix its name and path in the .in file
   3 - Run LAMMPS:
       $ mpirun -np xy <path>/lmp -in run_*.in -var seed $RANDOM
     (the seed is used to initialize the velocities, $RANDOM is a bash command)
  Your LAMMPS distribution needs to be compiled with the MOLECULE package.
  There is also a sample tk script that you can load in VMD (Visual Molecular
  Dynamics) to visualize your trajectory:
  $ vmd -e vmdscript.tk trajectoryfilename.lammpstrj


- printpotential: contains a C++ program that can be used to generate the
  tabulated potential that LAMMPS require. Build the program with
  $ g++ printPotential.cpp -o executable_name
  adapt the inputfile to your choice of parameters, then run with
  ./executable_name inputfile.txt output_directory
  (change executable_name and output_directory to your liking)


- startingstate_creators: contains Python scripts that can be used to generate
  different type of startingstates: cubic lattice, crystal-liquid layers,
  and single planes of different geometry plus cubic lattice that we use for
  the epitaxy runs. All of them have a help option, so run them with the -h
  flag to get explanations.


- postprocess: contains a C++ program to analize trajectories, requires CMake 3.5
  to build. There is a build script inside the directory; run it,
  the executable will appear in bld/lammpsIPCpostprocess.
  It computes the pair distribution function, the number of bonded neighbours,
  and an histogram of the orientations.
  Run the executable without any flag and the help will appear. The inputfile.txt
  that it asks is in the same directory, copy it and adapt it.
  IMPORTANT: at the moment it only supports symmetric IPCs.
