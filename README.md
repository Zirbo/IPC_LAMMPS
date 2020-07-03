IPC for LAMMPS
Implementation of Inverse Patchy Colloids for LAMMPS.

We provide sample .in files and startingstates for LAMMPS, and
a program to print a tabulate of the IPC potential, so that
anybody interested can start simulating IPCs in LAMMPS in few minutes :)

For any help, feel free to contact us!

LAST UPDATED: 3 July 2020

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
  Your LAMMPS distribution needs to be compiled with the MOLECULE package,
  and the RIGID package if you use rigid.
  There is also a sample tk script that you can load in VMD (Visual Molecular
  Dynamics) to visualize your trajectory:
  $ vmd -e vmdscript.tk trajectoryfilename.lammpstrj


- printpotential:
  contains programs that compute the potentials in LAMMPS format.
  Two scripts are given:
  -1- mapEpsilonsToCoefficients.sh
     requires as input (modify the file to supply the inputs!)
     the geometric parameters and the epsilons from the mapping (as the defined
     in the original Bianchi-Kahl-Likos paper) and prints out the potential
     in a lammpspot_name directory and a recap file inputfile_name.txt;
  -2- mapContactValuesToCoefficients.sh
     requires as input (modify the file to supply the inputs!)
     the geometric parameters and the desired contact values of the potential
     in the three EE, EP and PP orientations, and prints out the potential
     in a lammpspot_name directory, together with a file MC_inputfile_name.txt
     that recaps what you inputted and another file inputfile_name.txt (that
     follows the normalization of mapEpsilonsToCoefficients.sh)
  THE SCRIPTS ASSUME g++ AND python3 TO BE INSTALLED!


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
