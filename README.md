<h1>IPC for LAMMPS</h1>
<p>Implementation of Inverse Patchy Colloids for LAMMPS.</p>

<p>We provide sample .in files and startingstates for LAMMPS, and
a program to print a tabulate of the IPC potential, so that
anybody interested can start simulating IPCs in LAMMPS in few minutes :)</p>

<p>For any help, feel free to contact us!</p>

<p>LAST UPDATED: 6 Dec 2021</p>

<p>We have 4 directories.</p>

<ol>
  <li> make_potentials <br>
    contains programs that compute the potentials in LAMMPS format. <br>
    Two scripts are given:<br>
    <ul>
       <li>
       from_epsilon_to_potentials.sh <br>
       requires as input (modify the file to supply the inputs!)
       the geometric parameters and the epsilons from the mapping (as the defined
       in the original Bianchi-Kahl-Likos paper) and prints out the potential
       in a lammpspot_name directory and a recap file inputfile_name.txt;
       </li>
       <li>
       from_contact_values_to_potentials.sh <br>
       requires as input (modify the file to supply the inputs!)
       the geometric parameters and the desired contact values of the potential
       in the three EE, EP and PP orientations, and prints out the potential
       in a lammpspot_name directory, together with a file MC_inputfile_name.txt
       that recaps what you inputted and another file inputfile_name.txt (that
       follows the normalization of mapEpsilonsToCoefficients.sh)
       </li>
    </ul>
    IMPORTANT: you need a c+11-compatible version of g++ installed to build.
  </li>

  <li> lammps_inputfiles <br>
    each subdirectory of it contains a sample .infile and the
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
  </li>

  <li> startingstate_creators <br>
    contains Python3 scripts that can be used to generate FCC startingstates
    for ipcs and janus-ipcs. There's also one cubic lattice for FCC only.
    All of them have a help, so run them with the -h flag to get explanations.
  </li>

  <li> advanced <br>
    contains other startingstate generators, postprocess and visualization
    tools that we have been using in our research.
    You are welcome to use them, but they are not all very well documented,
    so use them at your own risk ;)
  </li>
</ol>
