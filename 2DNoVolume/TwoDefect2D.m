%Initial Condition for 2 uniaxial defects

function [Mesh,q,DoFMap] = TwoDefect2D(Sn,m1,m2,r1,r2,phi0,dphi,R,numtri)
   if (max(r1) > R) || (max(r2) > R)
       error('Defects must be located in domain');
   end
   
   %Create Mesh
   [Omega_Tri, Omega_Vertex] = bcc_triangle_mesh(numtri,numtri);
   %Go from -R to R
   Omega_Vertex = R*(2*Omega_Vertex - 1);
   
   Mesh = MeshTriangle(Omega_Tri,Omega_Vertex,'Omega');
   Bdy_Edges = Mesh.freeBoundary();
   Mesh = Mesh.Append_Subdomain('1D','bdy',Bdy_Edges);
   
   %Get DoFMap
   DoFMap = uint32(Mesh.ConnectivityList);
   
   %Create Initial Condition with linear profile for S
   S = Sn*ones(length(Mesh.Points),1);
   for ii = 1:length(Mesh.Points)
       if sqrt((Mesh.Points(ii,1) - r1(1))^2 + (Mesh.Points(ii,2) - r1(2))^2) < 0.1
           S(ii) = (10*Sn)*sqrt((Mesh.Points(ii,1) - r1(1))^2 + (Mesh.Points(ii,2) - r1(2))^2);
       elseif sqrt((Mesh.Points(ii,1) - r2(1))^2 + (Mesh.Points(ii,2) - r2(2))^2) < 0.1
           S(ii) = (10*Sn)*sqrt((Mesh.Points(ii,1) - r2(1))^2 + (Mesh.Points(ii,2) - r2(2))^2);
       end
   end
   
   %Create Initial Condition with Defects for phi
   r = norm(r1 - r2);
   phi = m1*atan2(Mesh.Points(:,2) - r1(2),Mesh.Points(:,1) - r1(1)) + m2*atan2(Mesh.Points(:,2) - r2(2),Mesh.Points(:,1) - r2(1)) + 0.5*dphi*(1 + (log((Mesh.Points(:,1) - r1(1)).^2 + (Mesh.Points(:,2) - r1(2)).^2) - log((Mesh.Points(:,1) - r2(1)).^2 + (Mesh.Points(:,2) - r2(2)).^2))./(log(r^2) - log(0.01))) + phi0;
   
   %Create Initial q
   q(:,1) = (1/sqrt(3)).*S.*(1 - 1.5*sin(phi).^2);
   q(:,2) = 0.5.*S.*sin(phi).^2;
   q(:,3) = 0.5*S.*sin(2*phi);
   q(:,4) = 0*Mesh.Points(:,1);
   q(:,5) = 0*Mesh.Points(:,1);
end