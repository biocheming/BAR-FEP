# BAR-FEP
BAR and FEP method for free energy calculation
**About BAR-FEP:**

This programe give the free energy calculation methods of BAR and FEP. The code is writen in Fortran. Here, we also give an example for calculating the solvation free energy of methane using BAR anf FEP. The input files for MD simulation are listed in the "in" directory.

**Code structure:**

The code is organized as follows:<br>
"./in": The input files for running MD simulation. In this example, five windows are set for FEP and BAR calculation.<br>
"./cal/Bar_Fep.f90": The main programe for BAR and FEP calculation.<br>
"./cal/*.dat": The MD output file from the step "prod" in every windows.<br>
"./cal/filename.dat": The name of the output dat file.<br>

**Use of this code**

set the windows and frames you used in Bar_Fep.f90 (windows=5, frames=10000 are set in this example) <br>
gfortran Bar_Fep.f90 <br>
./a.out
