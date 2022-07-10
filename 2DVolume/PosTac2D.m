%Initial condition for positive tactoid

function [Mesh,q,DoFMap] = PosTac2D(Sn,Rat,beta,R0,R,numtri)
    %Find tactoid lengths
    r2 = sqrt(0.5*(R0^2)*pi/(((((Rat^2 + 1)^2))/4)*atan2(2*Rat,Rat^2 - 1) - (Rat*(Rat^2 - 1))/2));
    r1 = Rat*r2;
    if r1 > R
        error('Tactoid does not fit in domain!');
    end
    r1n = beta*r1;
    Rcir = (r1^2 + r2^2)/(2*r2);
    Rcirn = (r1n^2 + r2^2)/(2*r2);
    y0 = (r1^2 - r2^2)/(2*r2);
    
    %Create Mesh
    [Omega_Tri,Omega_Vertex] = bcc_triangle_mesh(numtri,numtri);
    %Go from -R to R
    Omega_Vertex = R*(2*Omega_Vertex - 1);
    Mesh = MeshTriangle(Omega_Tri, Omega_Vertex,'Omega');
    Bdy_Edges = Mesh.freeBoundary();
    Mesh = Mesh.Append_Subdomain('1D','Bdy',Bdy_Edges);
    
    %Get DoFMap
    DoFMap = uint32(Mesh.ConnectivityList);
    
    %Create Initial Conditions on Mesh
    r = sqrt(Omega_Vertex(:,1).^2 + Omega_Vertex(:,2).^2);
    S = 0*Omega_Vertex(:,1);
    n = zeros(length(Mesh.Points(:,1)),2);
    for i=1:length(Omega_Vertex)
        if (Mesh.Points(i,1) > -r1) && (Mesh.Points(i,1) < r1) && (Mesh.Points(i,2) < sqrt(Rcir^2 - Mesh.Points(i,1)^2) - y0) && (Mesh.Points(i,2) > -sqrt(Rcir^2 - Mesh.Points(i,1)^2) + y0)
            S(i) = Sn;
            if Mesh.Points(i,2) > 0
                tfac = Mesh.Points(i,2)/(sqrt(Rcirn^2 - Mesh.Points(i,1)^2) - y0);
                n(i,:) = [sqrt(Rcirn^2 - Mesh.Points(i,1)^2),-tfac*Mesh.Points(i,1)]./sqrt(Rcirn^2 + (tfac^2 - 1)*Mesh.Points(i,1)^2);
            else
                tfac = Mesh.Points(i,2)/(-sqrt(Rcirn^2 - Mesh.Points(i,1)^2) + y0);
                n(i,:) = [sqrt(Rcirn^2 - Mesh.Points(i,1)^2),tfac*Mesh.Points(i,1)]./sqrt(Rcirn^2 + (tfac^2 - 1)*Mesh.Points(i,1)^2);
            end
        end
    end
    
    q(:,1) = (sqrt(3)/2)*S.*(n(:,1).^2 - (1/3));
    q(:,2) = 0.5*S.*n(:,2).^2;
    q(:,3) = S.*n(:,1).*n(:,2);
    q(:,4) = 0*Mesh.Points(:,1);
    q(:,5) = 0*Mesh.Points(:,1);
end