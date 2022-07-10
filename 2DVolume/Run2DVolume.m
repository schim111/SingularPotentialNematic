%-----Run2DVolume-----
%Setup file to run 2D Solver for Nematic with volume constraint
%Edit constants and initial condition here, then run this script to run
%code

%Setup Parallel pool within matlab for Lagrange Multipliers (not required)
%cl = parcluster('local');
%pool = cl.parpool(24);

%Add previous directory to path
addpath('..');

%Add Felicity to path
addpath(genpath('\FELICITY'));

%Add AGMG to path (not required)
%addpath(genpath('/AGMG'));

%Constants
%-----Material Params-----
alph = 3.4049; %Phase behavior of nematic
L = 1; %First Elastic Constant
LS = 0; %Special Elastic Constant, for adjusting surface tension only, not well tested
L2 = 0; %Second Elastic Constant
L3 = 0; %Third (cubic order) elastic constant
L4 = 0; %Fourth (cubic order) elastic constant

%-----Simulation Params-----
dt = 1e-1; %time-step
tol = 1e-6; %Energy stopping criteria (set to zero if not looking for equilibria)
numits = 50; %Maximum number of iterations
Dirichlet = 0; %Boundary conditions 0=Neumann, 1=Dirichlet
AG = 0; %Use AGMG? If 1 make sure AGMG is in the path
Constraint = 1; %If 0, the system will run without volume constraint

%-----Setup Lebedev Sphere-----
%degree: { 6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302, 
%        350, 434, 590, 770, 974, 1202, 1454, 1730, 2030, 2354, 2702, 3074, 
%        3470, 3890, 4334, 4802, 5294, 5810 } <--- Choose quadrature degree
Quad_Deg = 590;
LEB = getLebedevSphere(Quad_Deg);

%-----Find Sn-----
%Either manually enter or use function MSMin
%Sn = 0.4281;
[~,~,~,Sn,~] = MSMin(alph,LEB);

%Initial Conditions

%-----Load Mesh and q data----- <--May choose to load own initial q
%load('');

%-----Mesh Constants-----
R = 15; %Domain radium (Domain will be 2R x 2R square)
numtri = 150; %Number of triangles in the mesh (will use bcc mesh)

%-----Positive Tactoid-----

Rat = 1; %Aspect ratio of tactoid
beta = 1; %Ratio of virtual defects to length
R0 = 5; %Radius of circle of equivalent area to tactoid
[Mesh,q_0,DoFMap] = PosTac2D(Sn,Rat,beta,R0,R,numtri);

%-----Negative Tactoid-----
%{
m = 0.5; %Topological charge of boundary conditions
R0 = 5; %Radius of circle of equivalent area to tactoid
[Mesh,q_0,DoFMap] = NegTac2D(Sn,m,R0,R,numtri);
%}

%-----Run Solver-----
[S,Q,F] = SolverVolume2D(q_0,Sn,Mesh,DoFMap,LEB,alph,L,LS,L2,L3,L4,dt,tol,numits,Dirichlet,AG,Constraint);

%-----Setup Data-----
filename = ['','.mat']; %Data filename
save(['.\',filename],'Q','S','Mesh','F'); %Save to some directory














