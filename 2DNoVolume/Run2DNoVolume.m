%-----Run2DNoVolume-----
%Setup file to run 2D Solver for Nematic without volume constraint
%Edit constants and initial condition here, then run this script to run
%code

%Setup Parallel pool within matlab for Lagrange Multipliers (not required)
%cl = parcluster('local');
%pool = cl.parpool(24);

%Add previous directory to path
addpath('..');

%Add Felicity to path
addpath(genpath('/FELICITY'));

%Add AGMG to path (not required)
%addpath(genpath('/AGMG'));

%Constants
%-----Material Params-----
alph = 4; %Phase behavior of nematic
L = 1; %First Elastic Constant
L2 = 0; %Second Elastic Constant
L3 = 0; %Third (cubic order) elastic constant

%-----Simulation Params-----
dt = 1e-1; %time-step
tol = 1e-6; %Energy stopping criteria (set to zero if not looking for equilibria)
numits = 50; %Maximum number of iterations
Dirichlet = 0; %Boundary conditions 0=Neumann, 1=Dirichlet
AG = 0; %Use AGMG? If 1 make sure AGMG is in the path

%-----Setup Lebedev Sphere-----
%degree: { 6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302, 
%        350, 434, 590, 770, 974, 1202, 1454, 1730, 2030, 2354, 2702, 3074, 
%        3470, 3890, 4334, 4802, 5294, 5810 } <--- Choose quadrature degree
Quad_Deg = 590;
LEB = getLebedevSphere(Quad_Deg);

%-----Find Sn-----
%Either manually enter or use function MSMin
%Sn = 0.6751;
[~,~,~,Sn,~] = MSMin(alph,LEB);

%Initial Conditions
%-----Load Mesh and q data-----
%load('');

%-----Mesh Constants-----
R = 5; %Domain radius (Domain will be 2R x 2R square)
numtri = 150; %Number of triangles in the mesh (will use bcc mesh)

%-----Single Defect-----

m = 0.5; %Defect charge
r = [0,0]; %Defect Location
psi = 0; %Defect orientation
[Mesh,q_0,DoFMap] = SingleDefect2D(Sn,m,psi,r,R,numtri);

%-----Two Defects-----
%{
m1 = 0.5; %Defect 1 charge
m2 = -0.5; %Defect 2 charge
r1 = [-0.5*R,0]; %Location of defect 1
r2 = [0.5*R,0]; %Location of defect 2
phi0 = 0; %Overall Phase
dphi = 0; %Twistedness of defects
[Mesh,q_0,DoFMap = TwoDefect2D(Sn,m1,m2,r1,r2,phi0,dphi,R,numtri);
%}

%-----Field H-----
hphi = 0; %Azimuthal Angle of H
htheta = 0; %Colatitude of H
h = 0; %Magnitude of H
H = h*[cos(hphi)*sin(htheta)*ones(length(Mesh.Points)),sin(hphi)*sin(htheta)*ones(length(Mesh.Points)),cos(htheta)*ones(length(Mesh.Points))];

%-----Run Solver-----
[S,Q,F] = Solver2D(q_0,Mesh,DoFMap,LEB,alph,L,L2,L3,dt,tol,numits,H,Dirichlet,AG);

%-----Setup Data-----
filename = ['','.mat']; %Data filename
save(['.\',filename],'Q','S','Mesh','F'); %Save to some directory 





































