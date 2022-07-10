%-----Run3D-----
%Setup file to run 3D Solver for Nematic.
%Edit constants and initial condition here, then run this script to run the
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

%-----Mesh Constants-----
R = 5; %Domain radius (Domain will be 2R x 2R x 2R cube)
numpts = 150; %Number of pts along one edge of the mesh

%-----Single Disclination-----

p1 = [1,0,0]; %Initial vector defining n
OMEGA = [0,0,1]; %OMEGA of the line
T = [0,0,1]; %Tangent vector of disc
[Mesh,q_0,DoFMap] = SingleLine3D(Sn,p1,OMEGA,T,R,numpts);

%-----Single Zero Charge Loop-----
%{
p1 = [1,0,0]; %Initial vector defining n
OMEGA = [0,0,1]; %OMEGA of the loop
R0 = 2.5; %Initial radius of loop
[Mesh,q_0,DoFMap] = LoopDisc3D(Sn,p1,OMEGA,R0,R,numpts);
%}
%-----Two disclination lines-----
%{
p1 = [1,0,0]; %Initial vector defining n
OMEGA1 = [0,0,1]; %OMEGA of line 1
OMEGA2 = -[0,0,1]; %OMGEA of line 2
beta = 0; %Initial angle between lines
R0 = 5; %Initial separation of lines
[Mesh,q_0,DoFMap] = TwoLines3D(Sn,p1,OMEGA1,OMEGA2,beta,R0,R,numpts);
%}
%-----Saturn Ring-----
%{
R0 = 2.5; %Radius of Particle
[Mesh,q_0,DoFMap = SaturnRing(Sn,R0,R,numpts);
%}
%-----Velocity Field----- %Add static velocity field
gamm = 1; %Relaxation Velocity
xi = 1; %Tumbling parameter
V = zeros(length(Mesh.Points),3); %Velocity field

%Run Solver
[S,Q,F] = Solver3D(q_0,V,Mesh,DoFMap,LEB,alph,L,L2,L3,dt,tol,numits,Dirichlet,AG);

%-----Setup Data-----
filename = ['','.mat']; %Data filename
save(['.\',filename],'Q','S','Mesh','F'); %Save to some directory