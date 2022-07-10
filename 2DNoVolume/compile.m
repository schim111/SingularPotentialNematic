%Compile mex files
addpath(genpath('\FELICITY')); %<--- Set Path to FELICITY Directory
addpath('.\MatAssemblies');

Convert_Form_Definition_to_MEX(@MatAssem2D_MS1,{},'MS2DAssemble1');

Convert_Form_Definition_to_MEX(@MatAssem2D_MS2,{},'MS2DAssemble2');

Convert_Form_Definition_to_MEX(@MatAssem2D_MS3,{},'MS2DAssemble3');

Convert_Form_Definition_to_MEX(@MatAssem2D_MS4,{},'MS2DAssemble4');

Convert_Form_Definition_to_MEX(@MS_Energy_Assemble2D,{},'MS2DEnergy');
