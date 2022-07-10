%Calculates Maier-Saupe Energy for Biaxial Q
%Used to minimize the Maier-Saupe potential with respect to a Biaxial Q
function f = MaierSaupe(S,P,alph,LEB)
    %Get Lambda
    q = [(1/sqrt(3))*S,P,0,0,0];
    Lam = getLagMult(q,1e-8,LEB);
    %Find f
    Z = ZFunc(Lam,LEB);
    f = Lam(1)*(q(1) + sqrt(3)/6) + Lam(2)*(q(2) + 0.5) - log(Z) + log(4*pi) - 2*alph*(q(1)^2 + q(2)^2);
end