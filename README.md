# SingularPotentialNematic

Matlab code to run the Ball-Majumdar singular potential with anisotropic elasticity.
Included is code to run in two dimensions with/without a volume constraint and in three dimensions.
Also included are a few common initial conditions for configurations with defects or tactoids.

To run the code, the finite element, C++/Matlab, package FELICITY must be available 
(https://www.mathworks.com/matlabcentral/fileexchange/31141-felicity)

Included with the code is the file getLebedevSphere.m to get the quadrature points on a unit sphere.
(

Also optional is the AGMG solver for linear matrix equations
(http://agmg.eu/)

To use the code, first edit the compile.m file in the directory that is being used to include the 
path to FELCITY. Then run compile.m to create .mex files that will be used to assemble matrices.
If different equations are used, the matrix assembly files may be edited in the MatAssemblies
directories.

Then, edit the RunX.m file in the corresponding directory to include the appropriate material and 
simulation parameters and initial conditions. Finally, run the RunX.m file. The solutions will be 
saved where indicated in the RunX.m file.
