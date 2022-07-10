%Initial Condition for Negative Tactoid

function [Mesh,q,DoFMap] = NegTac2D(Sn,m,R0,R,numtri)
    if R0+15 > R
        error('Tactoid size too large! Must by < R - 15');
    end
    %Create Mesh
    [Omega_Tri,Omega_Vertex] = bcc_triangle_mesh(numtri,numtri);
    %Go from -R to R
    Omega_Vertex = R*(2*Omega_Vertex - 1);
    Mesh = MeshTriangle(Omega_Tri, Omega_Vertex,'Omega');
    Bdy_Edges = Mesh.freeBoundary();
    Mesh = Mesh.Append_Subdomain('1D','Bdy',Bdy_Edges);
    
    %Get DoF map
    DoFMap = uint32(Mesh.ConnectivityList);
    
    %Create initial conditions
    r = sqrt(Omega_Vertex(:,1).^2 + Omega_Vertex(:,2).^2);
    S = 0*Omega_Vertex(:,1);
    
    for ii=1:length(Omega_Vertex(:,1))
        if r(ii) > R0
            S(ii) = Sn;
        end
    end
    
    phi = m*atan2(Omega_Vertex(:,2),Omega_Vertex(:,1));
    
    for ii=1:length(Omega_Vertex(:,1))
        if r(ii) > R0 && r(ii) < R0 + 10
            phi(ii) = atan2(Omega_Vertex(ii,2),Omega_Vertex(ii,1)) + pi/2;
        end
    end
    
    q(:,1) = (1/sqrt(3))*(S - 1.5*S.*sin(phi).^2);
    q(:,2) = 0.5*S.*sin(phi).^2;
    q(:,3) = 0.5*S.*sin(2*phi);
    q(:,4) = 0*Omega_Vertex(:,1);
    q(:,5) = 0*Omega_Vertex(:,1);
end