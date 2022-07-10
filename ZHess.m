%{
function JFunc outputs the Jacobian matrix of the function Func defined in another file.
On input, JFunc takes a vector of lagrange multipliers, a vector of the order parameters,
and the Number of spaces for trapezoidal integration
%}
function M = ZHess(L,LEB)
    M = zeros(5,5);
    a = sqrt(3)/2; %Common multiplicitive factor
    A = (LEB.x).^2;
    B = (LEB.y).^2;
    C = (LEB.x).*(LEB.y);
    D = (LEB.x).*(LEB.z);
    E = (LEB.y).*(LEB.z);
    %Exponential factor
    g = a*L(1)*A + L(2)*(0.5*A + B) + L(3)*C + L(4)*D + L(5)*E;
    %integrands
    j11 = (a*A).*(a*A).*exp(g);
    j12 = (a*A).*(0.5*A + B).*exp(g);
    j13 = (a*A).*C.*exp(g);
    j14 = (a*A).*D.*exp(g);
    j15 = (a*A).*E.*exp(g);
    j22 = (0.5*A + B).*(0.5*A + B).*exp(g);
    j23 = (0.5*A + B).*C.*exp(g);
    j24 = (0.5*A + B).*D.*exp(g);
    j25 = (0.5*A + B).*E.*exp(g);
    j33 = (C).*C.*exp(g);
    j34 = (C).*D.*exp(g);
    j35 = (C).*E.*exp(g);
    j44 = (D).*D.*exp(g);
    j45 = (D).*E.*exp(g);
    j55 = (E).*E.*exp(g);
    %values of matrix
    M(1,1) = sum(j11.*LEB.w);
    M(1,2) = sum(j12.*LEB.w);
    M(1,3) = sum(j13.*LEB.w);
    M(1,4) = sum(j14.*LEB.w);
    M(1,5) = sum(j15.*LEB.w);
    M(2,2) = sum(j22.*LEB.w);
    M(2,3) = sum(j23.*LEB.w);
    M(2,4) = sum(j24.*LEB.w);
    M(2,5) = sum(j25.*LEB.w);
    M(3,3) = sum(j33.*LEB.w);
    M(3,4) = sum(j34.*LEB.w);
    M(3,5) = sum(j35.*LEB.w);
    M(4,4) = sum(j44.*LEB.w);
    M(4,5) = sum(j45.*LEB.w);
    M(5,5) = sum(j55.*LEB.w);
    M = M + M' - diag([M(1,1),M(2,2),M(3,3),M(4,4),M(5,5)]);
end
