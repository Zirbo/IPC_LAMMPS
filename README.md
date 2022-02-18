<h1>OSPC for LAMMPS</h1>
<p>Implementation of Overlapping Sphere Patchy Colloids for LAMMPS.</p>

<p>We provide all that is necessary to simulate OSPCs with LAMMPS:
sample .in files and startingstates, a program to print and tabulate the
OSPC potential, and generators of startingstates.
Anybody can start simulating IPCs in LAMMPS in few minutes.
For any help, feel free to contact us :)</p>

<p>We have 4 directories.</p>

<ol>
  <li> 1-make_potentials <br>
    contains programs that compute the potentials in LAMMPS format. <br>
    Three scripts are given:<br>
    <ul>
       <li>
       from_contact_values_to_potentials.sh <br>
       requires as input (modify the file to supply the inputs!)
       the geometric parameters and the desired contact values of the potential
       in the three EE, EP1 and P1P1 orientations
       (plus EP2, P1P2, P2P2 for asymmetric models)
       and prints out the potential in an output directory,
       named target_<modelname>_<symmetry>_contact.
       The directory contains a lammpspot dir that contains your potential.
       In addition there will be two directories,
       lammpspot_angular_plots and
       lammpspot_radial_plots,
       whose content can be plotted using gnuplot or xmgrace, so you can
       quickly visualize your potential. See the reference paper [2] for an
       explanation of the orientations.
       Finally there is also a file, inputfile.dat,
       that recaps your input values.
       </li>
       <li>
       from_epsilon_to_potentials.sh <br>
       requires as input (modify the file to supply the inputs!)
       the geometric parameters and the epsilons from the mapping (see the
       reference papers!) and prints out the potential in the output directory
       named target_<modelname>_<symmetry>_epsilons.
       The content is the same as for from_contact_values_to_potentials.sh.
       </li>
       <li>
       compute_ipc_geometry.py <br>
       in case you want to follow an ipc geometry, this python script helps
       you determine eccentricity and patch radius from patch amplitude and
       interaction range, or viceversa. check the inline manual with -h
       </li>
    </ul>
    IMPORTANT: you need a c+11-compatible version of g++ installed to be able to build.
  </li>

  <li> 2-lammps_inputfiles <br>
    each subdirectory of it contains a sample .in file and a suitable
    startingstate.
    <ul>
      <li>copy them to the directory where you want to run</li>
      <li>copy in that directory the potential that you want to use
        (see printpotential below!)
        and fix its name and path in the .in file</li>
      <li>run LAMMPS: <br>
         <code>$ mpirun -np <num_cpus> </lammps/exe/path> -in run_*.in -var seed $RANDOM </code> <br>
       where the seed is used to initialize the velocities, $RANDOM is a bash command <br>
    Your LAMMPS distribution needs to be compiled with the MOLECULE package.
    There is also a sample tk script that you can load in VMD (Visual Molecular
    Dynamics) to visualize your trajectory:<br>
    <code>$ vmd -e vmdscript.tk </path/to/trajectory.lammpstrj></code>
        </li>
    </ul>
  </li>

  <li> 3-startingstate_creators <br>
    contains Python3 scripts that can be used to generate FCC startingstates
    for janus and two-patch ospcs. There's also one cubic lattice for FCC only.
    All of them have a help, so run them with the -h flag to get explanations.
  </li>

  <li> advanced <br>
    contains other startingstate generators, postprocess and visualization
    tools that we have been using in our research.
    You are welcome to peep, but they are not all well documented,
    so use them at your own risk ;)
  </li>
</ol>

<p>We hope that you will have as much fun playing with OSPCs as we did :D</p>
<p>S.F. & E.B. </p>
