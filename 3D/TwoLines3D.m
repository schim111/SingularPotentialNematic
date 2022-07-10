%Initial Condition for two disclinations

function [Mesh,q,DoFMap] = TwoLines3D(Sn,p1,OMEGA1,OMEGA2,beta,R0,R,numpts)
    p1 = p1./norm(p1);
    OMEGA1 = OMEGA1./norm(OMEGA1);
    OMEGA2 = OMEGA2./norm(OMEGA2);
    if (abs(p1*OMEGA1') > 1e-5) || (abs(p1*OMEGA2') > 1e-5) 
        error('p1 must be orthogonal to OMEGA1 and OMEGA2!');
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
    for ii=1:M
        if sqrt((Mesh.Points(ii,1) + 0.5*R0)^2 + Mesh.Points(ii,2)^2) < 0.5
            S(ii) = 2*Sn*sqrt((Mesh.Points(ii,1) + 0.5*R0)^2 + Mesh.Points(ii,2)^2);
        elseif sqrt((-cos(beta)*Mesh.Points(ii,2) + sin(beta)*Mesh.Points(ii,3))^2 + (0.5*R0 - Mesh.Points(ii,1)^2)) < 0.5
            S(ii) = 2*Sn*sqrt((-cos(beta)*Mesh.Points(ii,2) + sin(beta)*Mesh.Points(ii,3))^2 + (0.5*R0 - Mesh.Points(ii,1)^2));
        end
    end
    
    %Setup director init
    p1 = [p1(1)*ones(M,1),p1(2)*ones(M,1),p1(3)*ones(M,1)];
    OMEGA1 = [OMEGA1(1)*ones(M,1),OMEGA1(2)*ones(M,2),OMEGA1(3)*ones(M,3)];
    OMEGA2 = [OMEGA2(1)*ones(M,1),OMEGA2(2)*ones(M,2),OMEGA2(3)*ones(M,3)];
    p21 = cross(OMEGA1,p1);
    p22 = cross(OMEGA2,p1);
    
    phi1 = atan2(Mesh.Points(:,2),Mesh.Points(:,1) + 0.5*R0);
    phi2 = atan2(-cos(beta)*Mesh.Points(:,2) + sin(beta)*Mesh.Points(:,3),0.5*R0 - Mesh.Points(:,1));
    
    n = (cos(0.5*phi1).*cos(0.5*phi2) - dot(OMEGA1,OMEGA2,2).*sin(0.5*phi1).*sin(0.5*phi2)).*p1 + sin(0.5*phi1).*cos(0.5*phi2).*p21 + sin(0.5*phi2).*cos(0.5*phi1).*p22;
    
    q(:,1) = (sqrt(3)/2)*S.*(n(:,1).*n(:,1) - (1/3));
    q(:,2) = 0.5*S.*(n(:,2).*n(:,2) - n(:,3).*n(:,3));
    q(:,3) = S.*n(:,1).*n(:,2);
    q(:,4) = S.*n(:,1).*n(:,3);
    q(:,5) = S.*n(:,2).*n(:,3);
end