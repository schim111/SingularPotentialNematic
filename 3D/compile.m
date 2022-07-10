%Compile mex files
addpath(genpath('.\FELICITY')); %<-- Add FELICITY to path
addpath('.\MatAssemblies');

Convert_Form_Definition_to_MEX(@VMatAssem3D_MS1,{},'VMS3DAssemble1');

Convert_Form_Definition_to_MEX(@VMatAssem3D_MS2,{},'VMS3DAssemble2');

Convert_Form_Definition_to_MEX(@VMatAssem3D_MS3,{},'VMS3DAssemble3');

Convert_Form_Definition_to_MEX(@VMatAssem3D_MS4,{},'VMS3DAssemble4');

Convert_Form_Definition_to_MEX(@MS_Energy_Assemble3D,{},'MS3DEnergy');
