%Initial Condition for a single disclination line

function [Mesh,q,DoFMap] = SingleLine3D(Sn,p1,OMEGA,T,R,numpts)
    T = T./norm(T);
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
    
    
    %Setup angle around disclination
    A = dot([1,0,0],T);
    B = dot([0,1,0],T);
    if A == 0
        a = 1;
        b = 0;
    else
        b = -1/sqrt((B^2)/(A^2) + 1);
        a = -(B/A)*b;
    end
    C = cross(T,[a,b,0]);
    phi = atan2(C(1)*Mesh.Points(:,1) + C(2)*Mesh.Points(:,2) + C(3)*Mesh.Points(:,3),a*Mesh.Points(:,1) + b*Mesh.Points(:,2));
    
    %Linearized core
    for ii=1:M
        if sqrt((a*Mesh.Points(ii,1) + b*Mesh.Points(ii,2))^2 + (C(1)*Mesh.Points(ii,1) + C(2)*Mesh.Points(ii,2) + C(3)*Mesh.Points(ii,3))^2) < 0.5
            S(ii) = 2*Sn*sqrt((a*Mesh.Points(ii,1) + b*Mesh.Points(ii,2))^2 + (C(1)*Mesh.Points(ii,1) + C(2)*Mesh.Points(ii,2) + C(3)*Mesh.Points(ii,3))^2);
        end
    end
    
    %Create Mesh scale vectors
    p1 = [p1(1)*ones(M,1),p1(2)*ones(M,1),p1(3)*ones(M,1)];
    OMEGA = [OMEGA(1)*ones(M,1),OMEGA(2)*ones(M,2),OMEGA(3)*ones(M,3)];
    p2 = cross(OMEGA,p1);
    n = cos(0.5*phi).*p1 + sin(0.5*phi).*p2;
    
    %Get q
    q(:,1) = (sqrt(3)/2)*S.*(n(:,1).*n(:,1) - (1/3));
    q(:,2) = 0.5*S.*(n(:,2).*n(:,2) - n(:,3).*n(:,3));
    q(:,3) = S.*n(:,1).*n(:,2);
    q(:,4) = S.*n(:,1).*n(:,3);
    q(:,5) = S.*n(:,2).*n(:,3);
end