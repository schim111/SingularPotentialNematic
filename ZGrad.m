%{
function Func outputs a vector of functions to be solved for in the overall
Maier-Saupe scheme. On input it takes the three values of the lagrange multipliers,
the three values of the order parameter field, and the number of points to use for the
trapezoid integration. LEB is a Lebedev set for quadrature on the unit
sphere
%}

function v = ZGrad(L,LEB)
    v = zeros(5,1);
    a = sqrt(3)/2; %Common multiplicitive factor 
    A = (LEB.x).^2;
    B = (LEB.y).^2;
    C = (LEB.x).*(LEB.y);
    D = (LEB.x).*(LEB.z);
    E = (LEB.y).*(LEB.z);
    %Exponential argument
    g = a*L(1)*A + L(2)*(0.5*A + B) + L(3)*C + L(4)*D + L(5)*E;
    %Integrands
    f1 = (a*A).*exp(g);
    f2 = (0.5*A + B).*exp(g);
    f3 = (C).*exp(g);
    f4 = (D).*exp(g);
    f5 = (E).*exp(g);
    %Trapezoid rule integrations
    v(1) = sum(f1.*LEB.w);
    v(2) = sum(f2.*LEB.w);
    v(3) = sum(f3.*LEB.w);
    v(4) = sum(f4.*LEB.w);
    v(5) = sum(f5.*LEB.w);
end
