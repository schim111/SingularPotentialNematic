function [S,Q,F_F] = SolverVolume2D(q_0,Sn,Mesh,DoFMap,LEB,alph,L,LS,L2,L3,L4,dt,tol,numits,Dirichlet,AG,Constraint)
    %Get DoFs For Triangle Mesh
    Soln_DoF = unique(Mesh.ConnectivityList(:));
    Num_Soln_DoF = length(Soln_DoF);
    Bdy_DoF = unique(Mesh.freeBoundary());
    q_F = zeros(length(Mesh.Points),5,numits+1);
    q_F(:,:,1) = q_0;
    F_F = zeros(numits+1,1);
    R = max(Mesh.Points(:,1));
    LSTemp = LS;
    L2Temp = L2;
    L3Temp = L3;
    L4Temp = L4;
    L2 = 0;
    L3 = 0;
    L4 = 0;
    
    %Get initial energy and Mu
    Lam = getLagMult(q_0,1e-9,LEB);
    Z = 0*Mesh.Points(:,1);
    for ii=1:length(Mesh.Points(:,1))
        Z(ii) = ZFunc(Lam(ii,:),LEB);
    end
    if ~Constraint
        mu = 0;
    else
        Mu = GetMu([],Mesh.Points,DoFmap,[],[],DoFmap,L,L2,L3,L4,LS,Lam,alph,q_0);
        V = (2*R)^2;
        mu = (Sn/(12*V))*Mu(1).MAT;
    end
    E = MS2DEnergyVol([],Mesh.Points,DoFMap,[],[],DoFMap,DoFMap,L,L2,L3,L4,LS,Lam,Z,alph,q_0);
    F0 = E(1).MAT;
    F_F(1) = F0;
    FEM1 = [];
    FEM2 = [];
    FEM3 = [];
    FEM4 = [];
    for its = 1:numits
        if mod(its,5) == 0
            disp(num2str(its));
        end
        
        if its==2 %Delay the application of Ls in case of instability
            LS = LSTemp;
            L2 = L2Temp;
            L3 = L3Temp;
            L4 = L4Temp;
        end
        
        %First half of iterations, get q_bar_1 and q_bar_2
        q_bar_1(:,1:2) = q_0(:,1:2);
        q_bar_2(:,1:3) = q_0(:,3:5);
        
        %Get Mats
        FEM1 = MS2DAssemble1Vol(FEM1,Mesh.Points,DoFMap,[],[],DoFMap,DoFMap,L,L2,L3,L4,LS,Sn,alph,dt,mu,q_0);
        LHSMass1 = FEM1(1).MAT;
        LHSStiff1 = FEM1(3).MAT;
        RHS1 = FEM1(2).MAT;
        %Newton sub iteration to solve linearized system for dq
        for newtits = 1:50 %usually converges in much less
            %Init dq
            dq = zeros(2*Num_Soln_DoF,1);
            
            %Get Boundary fixed nodes
            if Dirichlet %Dirichlet
                FixedNodes = [Bdy_DoF; Bdy_DoF + Num_Soln_DoF];
            else %Neumann
                FixedNodes = [];
            end
            FreeNodes = setdiff((1:1:length(dq))',FixedNodes);
            
            %Find Lam and dLam
            q = [q_bar_1,q_bar_2];
            [Lam, dLam] = getLagMult(q,1e-9,LEB);
            Lam = Lam(:,1:2);
            %Get Mats
            FEM2 = MS2DAssemble2Vol(FEM2,Mesh.Points,DoFMap,[],[],DoFMap,DoFMap,DoFMap,L,L2,L3,L4,LS,Lam,dLam,dt,q_bar_1,q_0);
            LHSMass2 = FEM2(1).MAT;
            RHS2 = FEM2(2).MAT;
            LHSStiff2 = FEM2(3).MAT;
            LHS = LHSMass1 + LHSStiff1 + LHSMass2 + LHSStiff2;
            RHS = RHS1 + RHS2;
            if AG
                dq(FreeNodes,1) = agmg(LHS(FreeNodes,FreeNodes),RHS(FreeNodes,1));
            else
                dq(FreeNodes,1) = LHS(FreeNodes,FreeNodes) \ RHS(FreeNodes,1);
            end
            
            %update q_bar
            q_bar_1(:,1) = q_bar_1(:,1) + dq(1:Num_Soln_DoF,1);
            q_bar_1(:,2) = q_bar_1(:,2) + dq(Num_Soln_DoF + 1:end,1);
            
            %Check inner Newton convergence
            err = max(abs(dq));
            if err < 1e-6
                break;
            end
            if newtits == 50
                error('Newton 1 did not converge!');
            end
        end
        
        %Update q_0
        q_0(:,1:2) = q_bar_1;
        
        %Second half of iterations
        %Get Mats
        FEM3 = MS2DAssemble3Vol(FEM3,Mesh.Points,DoFMap,[],[],DoFMap,DoFMap,L,L2,L3,L4,alph,dt,q_0);
        LHSMass1 = FEM3(1).MAT;
        LHSStiff1 = FEM3(1).MAT;
        RHS1 = FEM3(1).MAT;
        %Newton sub iteration to solve linearized system for dq
        for newtits = 1:50
            %Init dq
            dq = zeros(3*Num_Soln_DoF,1);
            
            %Get Boundary fixed nodes
            if Dirichlet %Dirichlet
                FixedNodes = [Bdy_DoF; Bdy_DoF + Num_Soln_DoF; Bdy_DoF + 2*Num_Soln_DoF];
            else %Neumann
                FixedNodes = [];
            end
            FreeNodes = setdiff((1:1:length(dq))',FixedNodes);
            
            %Get Lam and dLam
            q = [q_bar_1,q_bar_2];
            [Lam,~,dLam] = getLagMult(q,1e-9,LEB);
            Lam = Lam(:,3:5);
            
            %Get Mats
            FEM4 = MS2DAssemble4Vol(FEM4,Mesh.Points,DoFMap,[],[],DoFMap,DoFMap,DoFMap,L,L2,L3,L4,Lam,dLam,dt,q_bar_2,q_0);
            LHSMass2 = FEM4(1).MAT;
            RHS2 = FEM4(2).MAT;
            LHSStiff2 = FEM4(3).MAT;
            LHS = LHSMass1 + LHSStiff1 + LHSMass2 + LHSStiff2;
            RHS = RHS1 + RHS2;
            if AG
                dq(FreeNodes,1) = agmg(LHS(FreeNodes,FreeNodes),RHS(FreeNodes,1));
            else
                dq(FreeNodes,1) = LHS(FreeNodes,FreeNodes) \ RHS(FreeNodes,1);
            end
            
            %update q_bar
            q_bar_2(:,1) = q_bar_2(:,1) + dq(1:Num_Soln_DoF,1);
            q_bar_2(:,2) = q_bar_2(:,2) + dq(Num_Soln_DoF + 1:2*Num_Soln_DoF,1);
            q_bar_2(:,3) = q_bar_2(:,3) + dq(2*Num_Soln_DoF + 1:end,1);
            
            %Check convergence
            err = max(abs(dq));
            if err < 1e-6
                break;
            end
            if newtits == 50
                error('Newton 2 did not converge!');
            end
        end
        %update q_0
        q_0(:,3:5) = q_bar_2;
        q_F(:,:,its+1) = q_0;
        
        %Check energy change
        Lam = getLagMult(q_0,1e-9,LEB);
        Z = 0*Mesh.Points(:,1);
        for ii=1:length(Mesh.Points)
            Z(ii) = ZFunc(Lam(ii,:),LEB);
        end
        if Constraint
            Mu = GetMu2(Mu,Mesh.Points,DoFmap,[],[],DoFmap,L,L2,L3,L4,LS,Lam,alph,q_n);
            mu = (Sn/(12*V))*Mu(1).MAT;
        end
        E = MS2DEnergyVol([],Mesh.Points,DoFMap,[],[],DoFMap,DoFMap,L,L2,L3,L4,LS,Lam,Sn,Z,alph,q_0);
        F = E(1).MAT;
        F_F(its+1) = F;
        if abs(F - F0) < tol
            disp(['Energy convergence in ',num2str(its),' iterations']);
            q_F(:,:,its+2:end) = [];
            F_F(its+2:end) = [];
            break;
        end
        F0 = F;
    end
    
    %Get output S and Q
    s3 = 1/sqrt(3);
    S = zeros(length(Mesh.Points),size(q_F,3));
    Q = zeros(3,3,length(q_0),size(q_F,3));
    for ii=1:length(q_0)
        for kk=1:size(q_F,3)
            Q(:,:,ii,kk) = [2*s3*q_F(ii,1,kk),q_F(ii,3,kk),q_F(ii,4,kk);
                            q_F(ii,3,kk),-s3*q_F(ii,1,kk)+q_F(ii,2,kk),q_F(ii,5,kk);
                            q_F(ii,4,kk),q_F(ii,5,kk),-s3*q_F(ii,1,kk)-q_F(ii,2,kk)];
            S(ii,kk) = 1.5*max(eig(Q(:,:,ii,kk)));
        end
    end
end