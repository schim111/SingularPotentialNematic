%MSSolve2D script

%distcomp.feature('LocalUseMpiexec',false);
cl = parcluster('local');
pool = cl.parpool(24, 'IdleTimeout', Inf);



%Add Felicity to path
addpath(genpath('~/NematicIsoSoln/Matlab/FELICITY'));
addpath(genpath('~/NematicIsoSoln/Matlab/AGMG'));
rng(3);


%Define Constants: alph > 0, dt, L, L2, R, err, numits, numtri
m = 1;
alph = 4;
dt = 1e-1;
L = 1;
L2 = 0;
L3 = 0;
R = 5;
tol = 1e-8;
numits = 40;
numpts = 41;
beta = 0;


%get Lebedev sphere
% degree: { 6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302, 
%   350, 434, 590, 770, 974, 1202, 1454, 1730, 2030, 2354, 2702, 3074, 
%   3470, 3890, 4334, 4802, 5294, 5810 };
Quad_Deg = 1202;
LEB = getLebedevSphere(Quad_Deg);

%Get Sn from phase diagram (later can write method that calculates this)
Sn = 0.6751;
%phiD = -(5/27)*Sn*(L3/(2 + L2 + (2/3)*Sn*L3));
phiD = 0


%Create mesh with boundary subdomains
%Saturn Ring
%{
R2 = 1.5*7.5;
Cube_Dim = [-R, R];
Use_Newton = true;
TOL = 1e-12;
MG = Mesher3Dmex(Cube_Dim,numpts,Use_Newton,TOL);
LS = LS_Sphere();
LS.Param.cx = 0;
LS.Param.cy = 0;
LS.Param.cz = 0;
LS.Param.rad = 7.5;
LS.Param.sign = -1;
Interp_Handle = @(pt) LS.Interpolate(pt);
MG = MG.Get_Cut_Info(Interp_Handle);
[Omega_Tet,Omega_Vertex] = MG.run_mex(Interp_Handle);
clear MG;
clear LS;
clear Interp_Handle;
%}


[Omega_Tet,Omega_Vertex] = regular_tetrahedral_mesh(numpts,numpts,numpts);
%Go from -R to R
Omega_Vertex(:,1:2) = R*(2*Omega_Vertex(:,1:2) - 1);
Omega_Vertex(:,3) = R*(2*Omega_Vertex(:,3) - 1);

Mesh = MeshTetrahedron(Omega_Tet, Omega_Vertex,'Omega');
%Mesh = Mesh.Append_Subdomain('3D','NoHole',Omega_Tet);
%Mesh = Mesh.Output_Subdomain_Mesh('NoHole');
Bdy_Edges = Mesh.freeBoundary();
Mesh = Mesh.Append_Subdomain('2D','Bdy',Bdy_Edges);
Bdy_DoF = unique(Bdy_Edges(:));
Soln_DoF = unique(Mesh.ConnectivityList(:));
Num_Soln_DoF = length(Soln_DoF);
clear Omega_Tet;
clear Omega_Vertex;

%Find side boundary to fix sides
%{
FF_center = (1/3)*(Mesh.Points(Bdy_Edges(:,1),:) + Mesh.Points(Bdy_Edges(:,2),:) + Mesh.Points(Bdy_Edges(:,3),:));
Top_Mask = (FF_center(:,3) > 5-1e-5);
Bot_Mask = (FF_center(:,3) < -5+1e-5);
Top_Edges = Bdy_Edges(Top_Mask,:);
Bot_Edges = Bdy_Edges(Bot_Mask,:);
TopBot_Edges = [Top_Edges;Bot_Edges];
Side_Edges = setdiff(Bdy_Edges,TopBot_Edges,'stable','rows');
Mesh = Mesh.Append_Subdomain('2D','Side_Bdy',Side_Edges);
Side_DoF = unique(Side_Edges(:));
%}

%Get DoF map
DoFmap = uint32(Mesh.ConnectivityList);

%Create initial conditions q_0,q_n on mesh
%r = sqrt(Mesh.Points(:,1).^2 + Mesh.Points(:,2).^2 + Mesh.Points(:,3).^2);
%Domain_Centers = randi([1,Num_Soln_DoF],12,1);
%II = randi([1,1202],12,1);
S = Sn*ones(Num_Soln_DoF,1);
%{
for ii=Soln_DoF'
    if sqrt((Mesh.Points(ii,1) + 2)^2 + Mesh.Points(ii,2)^2) < 0.5
	S(ii) = 0;
    elseif sqrt((Mesh.Points(ii,1) - 2)^2 + (Mesh.Points(ii,3))^2) < 0.5
	S(ii) = 0;
    end
end
%}
%n = [zeros(Num_Soln_DoF,2),ones(Num_Soln_DoF,1)];
%sig = 1.25;
%for jj = 1:10
%    X = Mesh.Points(Domain_Centers(jj),:);
%    for ii=1:Num_Soln_DoF
%        r = norm(X - Mesh.Points(ii,:));
%        if r < 2*sig
	    %S(ii) = Sn*exp(-((r^2)/(2*sig^2)));
%	    S(ii) = Sn;
%	    n(ii,:) = [LEB.x(II(jj)),LEB.y(II(jj)),LEB.z(II(jj))];
%	end
%    end
%end



P = 0*Mesh.Points(:,1);
%Saturn Ring
%{
phi = atan2(Mesh.Points(:,2),Mesh.Points(:,1));
theta = atan2(sqrt(Mesh.Points(:,1).^2 + Mesh.Points(:,2).^2),Mesh.Points(:,3)) - 0.5*(atan2(sqrt(Mesh.Points(:,1).^2 + Mesh.Points(:,2).^2)-R2, Mesh.Points(:,3)) + pi/2);
%theta = 0*Mesh.Points(:,1);
%}

%Interacting Lines
%{
phi = 0*Mesh.Points(:,1);
theta = 0*Mesh.Points(:,1);
psi = 0*Mesh.Points(:,1);
for ii=Soln_DoF'
    if Mesh.Points(ii,1) < 0
	phi(ii) = -0.5*atan2(Mesh.Points(ii,2),Mesh.Points(ii,1)+2.5);
	theta(ii) = pi/2;
    elseif Mesh.Points(ii,1) > 0
	phi(ii) = 0.5*atan2(Mesh.Points(ii,2),Mesh.Points(ii,1)-2.5) + pi/2;
	theta(ii) = beta*(1 - sin(phi(ii)));
    else
	theta(ii) = pi/2;
    end
end
%}

%Saturn Ring Iniit
%{
for ii = Bdy_DoF'
    if r(ii) < 15
	theta(ii) = atan2(sqrt(Mesh.Points(ii,1)^2 + Mesh.Points(ii,2)^2),Mesh.Points(ii,3));
    else
	theta(ii) = 0;
    end
end

for ii = Soln_DoF'
    if r(ii) > 19
	theta(ii) = 0;
    end
end
%}
%
%R1 = 2.5;
%R2 = 4.5;
%
%r = sqrt(Mesh.Points(:,1).^2 + Mesh.Points(:,2).^2);
%phi = atan2(Mesh.Points(:,2),Mesh.Points(:,1));
n = zeros(Num_Soln_DoF,3);
%Skyrmions
%rhat = [cos(phi),sin(phi),zeros(Num_Soln_DoF,1)];
%zhat = [0,0,1];
%for ii=1:Num_Soln_DoF
%if r(ii) < 3
%    n(ii,:) = sin(pi*r(ii,:)/3).*rhat(ii,:) + cos(pi*r(ii,:)/3).*zhat;
%else
%    n(ii,:) = zhat;
%end
%end


p1 = [ones(Num_Soln_DoF,1),zeros(Num_Soln_DoF,1),zeros(Num_Soln_DoF,1)];
OM1 = -[zeros(Num_Soln_DoF,1),zeros(Num_Soln_DoF,1),ones(Num_Soln_DoF,1)];
p21 = cross(OM1,p1);
OM2 = [zeros(Num_Soln_DoF,1),sin(beta)*ones(Num_Soln_DoF,1),cos(beta)*ones(Num_Soln_DoF,1)];
p22 = cross(OM2,p1);
%phi = atan2(Mesh.Points(:,2) - Mesh.Points(:,3),Mesh.Points(:,1));
phi2 = atan2(-Mesh.Points(:,2),2 - Mesh.Points(:,1));
phi1 = atan2(Mesh.Points(:,2),Mesh.Points(:,1) + 2);
m = cos(0.5*phi1).*p1 + sin(0.5*phi1).*p21;
n = cos(0.5*phi2).*m + sin(0.5*phi2).*cross(OM2,m) + (1 - cos(0.5*phi2)).*dot(OM2,m).*OM2;
%n = cos(phi1 + phi2).*p1 + sin(phi1 + phi2).*p21;
%n = cos(0.5*atan2(Mesh.Points(:,2), - 2.5)).*p1 + sin(0.5*atan2(Mesh.Points(:,3),sqrt(Mesh.Points(:,1).^2 + Mesh.Points(:,2).^2) - 2.5)).*p2;
%{
n = zeros(Num_Soln_DoF,3);
for ii=Soln_DoF'
    if Mesh.Points(ii,1) < 1
	n(ii,:) = cos(0.5*atan2(Mesh.Points(ii,2),Mesh.Points(ii,1) + 2.5))*p1(ii,:) + sin(0.5*atan2(Mesh.Points(ii,2),Mesh.Points(ii,1) + 2.5))*p21(ii,:);
    elseif Mesh.Points(ii,1) > 1
	n(ii,:) = cos(0.5*atan2(Mesh.Points(ii,2),Mesh.Points(ii,1) - 2.5) + pi/2)*p1(ii,:) + sin(0.5*atan2(Mesh.Points(ii,2),Mesh.Points(ii,1) - 2.5) + pi/2)*p22(ii,:);
    else
	n(ii,:) = [1,0,0];
    end
end
%}

%theta = (pi)*rand(Num_Soln_DoF,1);
%phi = 2*pi*rand(Num_Soln_DoF,1);
%psi = 0*Mesh.Points(:,1);

%{
q_n(:,1) = (sqrt(3)/2)*(S.*((cos(phi).^2).*(sin(theta).^2) - (1/3)) + P.*(((sin(phi).^2) - (cos(phi).^2).*(cos(theta).^2)).*cos(2*psi) + sin(2*phi).*cos(theta).*sin(2*psi)));
q_n(:,2) = 0.125*(-2*(cos(phi).^2).*(S - 3*P.*cos(2*psi)) + cos(2*theta).*(-3 + cos(2*phi)).*(S + P.*cos(2*psi)) - 4*P.*cos(theta).*sin(2*phi).*sin(2*psi));
q_n(:,3) = S.*cos(phi).*sin(phi).*sin(theta).^2 - P.*(0.25*sin(2*phi).*(3 + cos(2*theta)).*cos(2*psi) + cos(2*phi).*cos(theta).*sin(2*psi));
q_n(:,4) = sin(theta).*((S + P.*cos(2*psi)).*cos(phi).*cos(theta) - P.*sin(phi).*cos(2*psi));
q_n(:,5) = sin(theta).*((S + P.*cos(2*psi)).*sin(phi).*cos(theta) + P.*cos(phi).*sin(2*psi));
%}
%
q_n(:,1) = (sqrt(3)/2)*S.*(n(:,1).*n(:,1) - (1/3));
q_n(:,2) = 0.5*S.*(n(:,2).*n(:,2) - n(:,3).*n(:,3));
q_n(:,3) = S.*n(:,1).*n(:,2);
q_n(:,4) = S.*n(:,1).*n(:,3);
q_n(:,5) = S.*n(:,2).*n(:,3);
%
q_0 = q_n;
q_F = zeros(Num_Soln_DoF,5,numits+1);
F_F = zeros(numits+1,1);
q_F(:,:,1) = q_n;

%Get unique DoFs
%Soln_DoF = unique(Mesh.ConnectivityList(:));
%Num_Soln_DoF = length(Soln_DoF);


%Get initial energy
Lam = getLagMult(q_n,1e-9,LEB);
Z = 0*Mesh.Points(:,1);
for j=1:Num_Soln_DoF
    Z(j) = ZFunc(Lam(j,:),LEB);
end
E = MS3DEnergy([],Mesh.Points,DoFmap,[],[],DoFmap,DoFmap,L,L2,L3,Lam,Z,alph,q_n);
F0 = E(1).MAT;
F_F(1) = F0;
FEM1 = [];
FEM2 = [];
FEM3 = [];
FEM4 = [];
%Run Time Loop
disp('Start Iterations');
for its = 1:numits
    disp(num2str(its));
    
    
    %Create q_bar_1 and q_bar_2 on mesh
    q_bar_1(:,1:2) = q_n(:,1:2);
    q_bar_2(:,1:3) = q_n(:,3:5);
    %First half time loop
    FEM1 = MS3DAssemble1(FEM1,Mesh.Points,DoFmap,[],[],DoFmap,DoFmap,L,L2,L3,alph,dt,q_n);
    LHSMass1 = FEM1(1).MAT;
    LHSStiff1 = FEM1(3).MAT;
    RHS1 = FEM1(2).MAT;
    %Newton sub iteration to solve linearized system for dq
    for newtits = 1:50 %placeholder number here
        dq = zeros(2*Num_Soln_DoF,1);
        %Get FreeNodes
        %FixedNodes = [Bdy_DoF; Bdy_DoF + Num_Soln_DoF];
	FixedNodes = [];
	%FixedNodes = [Side_DoF; Side_DoF + Num_Soln_DoF];
        FreeNodes = setdiff((1:1:length(dq))',FixedNodes);
        
        %Find Lam and dLam 
        q = [q_bar_1,q_bar_2];
        [Lam, dLam] = getLagMult(q,1e-9,LEB);
        Lam = Lam(:,1:2);
        FEM2 = MS3DAssemble2(FEM2,Mesh.Points,DoFmap,[],[],DoFmap,DoFmap,DoFmap,L,L2,L3,Lam,dLam,dt,q_bar_1,q_n);
        LHSMass2 = FEM2(1).MAT;
        RHS2 = FEM2(2).MAT;
        LHSStiff2 = FEM2(3).MAT;
        LHS = LHSMass1 + LHSStiff1 + LHSMass2 + LHSStiff2;
        RHS = RHS1 + RHS2;
        dq(FreeNodes,1) = agmg(LHS(FreeNodes,FreeNodes),RHS(FreeNodes,1));
        if norm(LHS(FreeNodes,FreeNodes)*dq(FreeNodes,1) - RHS(FreeNodes,1)) > 1e-5
	    disp(['First agmg nonsolve: |Ax-b| = ',num2str(norm(LHS(FreeNodes,FreeNodes)*dq(FreeNodes,1) - RHS(FreeNodes,1)))]);
	end
        %q_bar -> q_bar + dq
        q_bar_1(:,1) = q_bar_1(:,1) + dq(1:Num_Soln_DoF,1);
        q_bar_1(:,2) = q_bar_1(:,2) + dq(Num_Soln_DoF + 1:end,1);
        %Check convergence
        err = max(abs(dq));
        if err < 1e-5
            break;
        end
        if newtits == 50
            disp('Newton did not converge!');
        end
    end
    
    %q_n(1), q_n(2) -> q_bar
    q_n(:,1:2) = q_bar_1;
    
    %Second half time loop
    %Generate matrices 3
    FEM3 = MS3DAssemble3(FEM3,Mesh.Points,DoFmap,[],[],DoFmap,DoFmap,L,L2,L3,alph,dt,q_n);
    LHSMass1 = FEM3(1).MAT;
    LHSStiff1 = FEM3(3).MAT;
    RHS1 = FEM3(2).MAT;
    %Newton sub iteration solve linearized system for dq
    for newtits = 1:50 %placeholder number
        dq = zeros(3*Num_Soln_DoF,1);
        %Get free nodes
        %FixedNodes = [Bdy_DoF; Bdy_DoF + Num_Soln_DoF; Bdy_DoF + 2*Num_Soln_DoF];
	FixedNodes = [];
	%FixedNodes = [Side_DoF;Side_DoF + Num_Soln_DoF;Side_DoF + 2*Num_Soln_DoF];
        FreeNodes = setdiff((1:1:length(dq))',FixedNodes);
        
        %Get Lam and dLam
        q = [q_bar_1,q_bar_2];
        [Lam,~,dLam] = getLagMult(q,1e-9,LEB);
        Lam = Lam(:,3:5);
        
        %Generate matrices 4
        FEM4 = MS3DAssemble4(FEM4,Mesh.Points,DoFmap,[],[],DoFmap,DoFmap,DoFmap,L,L2,L3,Lam,dLam,dt,q_bar_2,q_n);
        LHSMass2 = FEM4(1).MAT;
        RHS2 = FEM4(2).MAT;
	LHSStiff2 = FEM4(3).MAT;
        LHS = LHSMass1 + LHSStiff1 + LHSMass2 + LHSStiff2;
        RHS = RHS1 + RHS2;
        dq(FreeNodes,1) = agmg(LHS(FreeNodes,FreeNodes),RHS(FreeNodes,1));
	if norm(LHS(FreeNodes,FreeNodes)*dq(FreeNodes,1) - RHS(FreeNodes,1)) > 1e-5
            disp(['Second agmg nonsolve: |Ax-b| = ',num2str(norm(LHS(FreeNodes,FreeNodes)*dq(FreeNodes,1) - RHS(FreeNodes,1)))]);
	end
	
        q_bar_2(:,1) = q_bar_2(:,1) + dq(1:Num_Soln_DoF,1);
        q_bar_2(:,2) = q_bar_2(:,2) + dq(Num_Soln_DoF + 1:2*Num_Soln_DoF,1);
        q_bar_2(:,3) = q_bar_2(:,3) + dq(2*Num_Soln_DoF + 1:end,1);
        
        %Check convergence
        err = max(abs(dq));
        if err < 1e-5
            break;
        end
        if newtits == 50
            disp('Newton did not converge!');
        end
    end
    q_n(:,3:5) = q_bar_2;
    q_F(:,:,its+1) = q_n;
    %check energy change, if small stop
    Lam = getLagMult(q_n,1e-9,LEB);
    Z = 0*Mesh.Points(:,1);
    for j=1:Num_Soln_DoF
        Z(j) = ZFunc(Lam(j,:),LEB);
    end
    E = MS3DEnergy(E,Mesh.Points,DoFmap,[],[],DoFmap,DoFmap,L,L2,L3,Lam,Z,alph,q_n);
    F = E(1).MAT;
    F_F(its+1) = F;
    if abs(F - F0) < tol*(abs(F) + abs(F0) + eps)
        disp(['Convergence in ',num2str(its),' iterations']);
	q_F(:,:,its+2:end) = [];
	F_F(its+2:end) = [];
        break;
    end
    F0 = F;
end

%Create Q
S = zeros(Num_Soln_DoF,size(q_F,3));
Q = zeros(3,3,length(q_n),size(q_F,3));
for ii=1:length(q_n)
    for k = 1:size(q_F,3);
        Q(:,:,ii,k) = [(2/sqrt(3))*q_F(ii,1,k),q_F(ii,3,k),q_F(ii,4,k);
                        q_F(ii,3,k),-(1/sqrt(3))*q_F(ii,1,k)+q_F(ii,2,k),q_F(ii,5,k);
                        q_F(ii,4,k),q_F(ii,5,k),-(1/sqrt(3))*q_F(ii,1,k)-q_F(ii,2,k)];
	S(ii,k) = 1.5*max(eig(Q(:,:,ii,k)));
    end
end

filename = ['InteractingLinesC10_L2=',num2str(L2/L),'_L3=',num2str(L3/L),'.mat'];
save(['Solns/',filename],'Q','S','Mesh','F_F');
    
    
    
    
        

