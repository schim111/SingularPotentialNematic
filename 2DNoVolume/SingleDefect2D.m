%Initial Condition for a uniaxial defect

function [Mesh,q,DoFMap] = SingleDefect2D(Sn,m,psi,r,R,numtri)
    if max(r) > R
       error('max(r) > R, Defect must be located in domain'); 
    end

    %Create Mesh
    [Omega_Tri,Omega_Vertex] = bcc_triangle_mesh(numtri,numtri);
    %Go from -R to R
    Omega_Vertex = R*(2*Omega_Vertex - 1);
    
    Mesh = MeshTriangle(Omega_Tri,Omega_Vertex,'Omega');
    Bdy_Edges = Mesh.freeBoundary();
    Mesh = Mesh.Append_Subdomain('1D','Bdy',Bdy_Edges);
    
    %Get DoFMap
    DoFMap = uint32(Mesh.ConnectivityList);
    
    %Create Initial Condition with linear profile for S
    S = Sn*ones(length(Mesh.Points),1);
    for ii=1:length(Mesh.Points)
        if sqrt((Mesh.Points(ii,1) - r(1))^2 + (Mesh.Points(ii,2) - r(2))^2) < 0.1
            S(ii) = (10*Sn)*sqrt((Mesh.Points(ii,1) - r(1))^2 + (Mesh.Points(ii,2) - r(2))^2);
        end
    end
    
    %Create Initial Condition with Defected Profile for phi
    phi = m*atan2(Mesh.Points(:,2) - r(2),Mesh.Points(:,1) - r(1)) + abs(m)*psi;
    
    %Create Initial q
    q(:,1) = (1/sqrt(3)).*S.*(1 - 1.5*sin(phi).^2);
    q(:,2) = 0.5.*S.*sin(phi).^2;
    q(:,3) = 0.5*S.*sin(2*phi);
    q(:,4) = 0*Mesh.Points(:,1);
    q(:,5) = 0*Mesh.Points(:,1);
end