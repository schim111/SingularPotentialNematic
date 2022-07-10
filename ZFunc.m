%{
function Func outputs a vector of functions to be solved for in the overall
Maier-Saupe scheme. On input it takes the three values of the lagrange multipliers,
the three values of the order parameter field, and the number of points to use for the
trapezoid integration.
%}

function v = ZFunc(L,LEB)
    A = (LEB.x).^2;
    B = (LEB.y).^2;
    C = (LEB.x).*(LEB.y);
    D = (LEB.x).*(LEB.z);
    E = (LEB.y).*(LEB.z);
    %Exponential argument
    g = (sqrt(3)/2)*L(1)*A + L(2)*(0.5*A + B) + L(3)*C + L(4)*D + L(5)*E;
    %Integrands
    f = exp(g);
    %Trapezoid rule integrations
    v = sum(f.*LEB.w);
end
