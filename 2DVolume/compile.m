%Compile mex files
addpath(genpath('\FELICITY')); %<---- Add Felicity to path
addpath('.\MatAssemblies');

Convert_Form_Definition_to_MEX(@MatAssem2D_MS1_Vol,{},'MS2DAssemble1Vol');

Convert_Form_Definition_to_MEX(@MatAssem2D_MS2_Vol,{},'MS2DAssemble2Vol');

Convert_Form_Definition_to_MEX(@MatAssem2D_MS3_Vol,{},'MS2DAssemble3Vol');

Convert_Form_Definition_to_MEX(@MatAssem2D_MS4_Vol,{},'MS2DAssemble4Vol');

Convert_Form_Definition_to_MEX(@MS_Energy_Assemble2D_Vol,{},'MS2DEnergyVol');

Convert_Form_Definition_to_MEX(@Get_Mu_Assemble,{},'GetMu');
