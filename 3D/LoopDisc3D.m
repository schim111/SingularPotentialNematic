%Initial Condition for a single zero charge disclination loop

function [Mesh,q,DoFMap] = LoopDisc3D(Sn,p1,OMEGA,R0,R,numpts)
    p1 = p1./norm(p1);
    OMEGA = OMEGA./norm(OMEGA);
    if abs(p1*OMEGA') > 1e-5
        error('p1 must be orthogonal to OMEGA!');
    end
    %Create Mesh
    [Omega_Tet,Omega_Vertex] = regular_tetrahedral_mesh(numpts,numpts,numpts);
    %Go from -R to R
    Omega_Vertex = R*(2*Omega_Vertex - 1);
    Mesh = MeshTetrahedron(Omega_Tet, Omega_Vertex,'Omega');
    
    %Boudnary
    Bdy_Edges = Mesh.freeBoundary();
    Mesh = Mesh.Append_Subdomain('2D','Bdy',Bdy_Edges);
    
    %Get length of Mesh
    M = length(Mesh.Points);
    
    %Get DoFMap
    DoFMap = unint32(Mesh.ConnectivityList);
    
    %Setup initial conditions with linearized core
    S = Sn*ones(M,1);
    
    %Setup director
    r = sqrt(Mesh.Points(:,1).^2 + Mesh.Points(:,2).^2);
    p1 = [p1(1)*ones(M,1),p1(2)*ones(M,1),p1(3)*ones(M,1)];
    OMEGA = [OMEGA(1)*ones(M,1),OMEGA(2)*ones(M,2),OMEGA(3)*ones(M,3)];
    p2 = cross(OMEGA,p1);
    phi1 = atan2(Mesh.Points(:,3),R0 - r);
    phi2 = atan2(-Mesh.Points(:,3),R0 + r);
    
    %Get n
    n = cos(0.5*(phi1 + phi2)).*p1 + sin(0.5*(phi1 + phi2)).*p2;
    
    %Linearize core
    for ii=1:M
        if sqrt((r(ii) - R0)^2 + Mesh.Points(ii,3)^2) < 0.5
            S(ii) = 2*Sn*sqrt((r(ii) - R0)^2 + Mesh.Points(ii,3)^2);
        end
    end
    
    q(:,1) = (sqrt(3)/2)*S.*(n(:,1).*n(:,1) - (1/3));
    q(:,2) = 0.5*S.*(n(:,2).*n(:,2) - n(:,3).*n(:,3));
    q(:,3) = S.*n(:,1).*n(:,2);
    q(:,4) = S.*n(:,1).*n(:,3);
    q(:,5) = S.*n(:,2).*n(:,3);
end