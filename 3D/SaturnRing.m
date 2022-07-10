%Initial Condition for a Saturn ring

function [Mesh,q,DoFMap] = SaturnRing(Sn,R0,R,numpts)
    %Create Mesh Using Felicity's Mesh generator
    R2 = 1.5*R0;
    Cube_Dim = [-R, R];
    Use_Newton = true;
    TOL = 1e-12;
    MG = Mesher3Dmex(Cube_Dim,numpts,Use_Newton,TOL);
    LS = LS_Sphere();
    LS.Param.cx = 0;
    LS.Param.cy = 0;
    LS.Param.cz = 0;
    LS.Param.rad = R0;
    LS.Param.sign = -1;
    Interp_Handle = @(pt) LS.Interpolate(pt);
    MG = MG.Get_Cut_Info(Interp_Handle);
    [Omega_Tet,Omega_Vertex] = MG.run_mex(Interp_Handle);
    
    Mesh = MeshTetrahedron(Omega_Tet, Omega_Vertex,'Omega');
    Mesh = Mesh.Remove_Unused_Vertices;
    Bdy_Edges = Mesh.freeBoundary();
    Mesh = Mesh.Append_Subdomain('2D','Bdy',Bdy_Edges);
    Bdy_DoF = unique(Bdy_Edges(:));
    
    DoFMap = uint32(Mesh.ConnectivityList);
    
    M = length(Mesh.Points);
    
    %Create Initial conditions with linearized core
    S = Sn*ones(M,1);
    
    r = sqrt(Mesh.Points(:,1).^2 + Mesh.Points(:,2).^2);
    
    phi1 = atan2(Mesh.Points(:,3),R2 - r);
    phi2 = atan2(-Mesh.Points(:,3),R2 + r);
    
    %Get n
    p1 = [Mesh.Points(:,1),Mesh.Points(:,2),zeros(M,1)]./r;
    p2 = [zeros(M,1),zeros(M,1),-ones(M,1)];
    n = sin(0.5*(phi1 + phi2)).*p1 + sin(0.5*(phi1 + phi2)).*p2;
    for ii=Bdy_DoF'
        if sqrt(Mesh.Points(ii,1)^2 + Mesh.Points(ii,2)^2 + Mesh.Points(ii,3)^2) < R2
            n(ii,:) = [Mesh.Points(ii,1),Mesh.Points(ii,2),Mesh.Points(ii,3)]./sqrt(Mesh.Points(ii,1)^2 + Mesh.Points(ii,2)^2 + Mesh.Points(ii,3)^2);
        else
            n(ii,:) = [0,0,1];
        end
    end
    
    for ii=1:M
        if sqrt((r(ii) - R2)^2 + Mesh.Points(ii,3)^2) < 0.5
            S(ii) = 2*Sn*sqrt((r(ii) - R2)^2 + Mesh.Points(ii,3)^2);
        end
    end
    
    q(:,1) = (sqrt(3)/2)*S.*(n(:,1).*n(:,1) - (1/3));
    q(:,2) = 0.5*S.*(n(:,2).*n(:,2) - n(:,3).*n(:,3));
    q(:,3) = S.*n(:,1).*n(:,2);
    q(:,4) = S.*n(:,1).*n(:,3);
    q(:,5) = S.*n(:,2).*n(:,3);
end
