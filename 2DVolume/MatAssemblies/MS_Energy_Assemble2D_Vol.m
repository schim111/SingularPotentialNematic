function MATS = MS_Energy_Assemble2D()

    %Define Domain
    Omega = Domain('triangle');
    
    %Define finite element spaces
    Vector_P1 = Element(Omega,lagrange_deg1_dim2,5);
    Scalar_P1 = Element(Omega,lagrange_deg1_dim2,1);
    
    %Functions
    q = Coef(Vector_P1);
    Lam = Coef(Vector_P1);
    Z = Coef(Scalar_P1);
    alph = Constant(Omega);
    L = Constant(Omega);
    LS = Constant(Omega);
    L2 = Constant(Omega);
    L3 = Constant(Omega);
    L4 = Constant(Omega);
    Sn = Constant(Omega);
    C1 = 4*pi;
    
    %Define eval matrices
    Fb = Real(1,1);
    Fe = Real(1,1);
    Ib = Integral(Omega,Lam.val(1).*(q.val(1) + (sqrt(3)/6)) + Lam.val(2).*(q.val(2) + 0.5) + Lam.val(3).*q.val(3) + Lam.val(4).*q.val(4) + Lam.val(5).*q.val(5) + log(C1) - log(Z.val) - 2*alph.val.*(q.val(1).^2 + q.val(2).^2 + q.val(3).^2 + q.val(4).^2 + q.val(5).^2));
    Ie = Integral(Omega,2*L.val*((q.grad(1)' * q.grad(1)) + (q.grad(2)' * q.grad(2)) + (q.grad(3)' * q.grad(3)) + (q.grad(4)' * q.grad(4)) + (q.grad(5)' * q.grad(5))) + LS.val*(3*(q.grad(1)' * q.grad(1)) + 9*(q.grad(2)' * q.grad(2)) + 6*sqrt(3)*(q.grad(1)' * q.grad(2))) + L2.val*((q.grad(3,2) + (2/sqrt(3))*q.grad(1,1))^2 + ((-1/sqrt(3))*q.grad(1,2) + q.grad(2,2) + q.grad(3,1))^2 + (q.grad(4,1) + q.grad(5,2))^2) + 2*L3.val*(-(1/sqrt(3))*q.val(1)*(q.grad(1,2)^2 + q.grad(2,2)^2 + q.grad(3,2)^2 + q.grad(4,2)^2 + q.grad(5,2)^2 - 2*(q.grad(1,1)^2 + q.grad(2,1)^2 + q.grad(3,1)^2 + q.grad(4,1)^2 + q.grad(5,1)^2)) + q.val(2)*(q.grad(1,2)^2 + q.grad(2,2)^2 + q.grad(3,2)^2 + q.grad(4,2)^2 + q.grad(5,2)^2) + 2*q.val(3)*(q.grad(1,1)*q.grad(1,2) + q.grad(2,1)*q.grad(2,2) + q.grad(3,1)*q.grad(3,2) + q.grad(4,1)*q.grad(4,2) + q.grad(5,1)*q.grad(5,2))) + (L4.val/3)*(-2*sqrt(3)*q.val(3)*q.grad(1,2)*q.grad(3,2) + 6*q.val(3)*q.grad(2,2)*q.grad(3,2) - 2*sqrt(3)*q.val(4)*q.grad(1,2)*q.grad(4,2) + 6*q.val(4)*q.grad(2,2)*q.grad(4,2) - 2*sqrt(3)*q.val(5)*q.grad(1,2)*q.grad(5,2) + 6*q.val(5)*q.grad(2,2)*q.grad(5,2) + 12*q.val(3)*q.grad(3,2)*q.grad(3,1) + 6*q.val(4)*q.grad(4,2)*q.grad(3,1) + 6*q.val(5)*q.grad(5,2)*q.grad(3,1) + 4*sqrt(3)*q.val(3)*q.grad(1,1)*q.grad(3,1) + q.val(1)*(-2*sqrt(3)*q.grad(1,2)^2 + 2*q.grad(1,1)*(3*q.grad(3,2) + 2*sqrt(3)*q.grad(1,1)) + 6*q.grad(1,2)*(q.grad(2,2) + q.grad(3,1))) + q.val(2)*(2*q.grad(2,1)*(3*q.grad(3,2) + 2*sqrt(3)*q.grad(1,1)) + 2*q.grad(2,2)*(-sqrt(3)*q.grad(1,2) + 3*(q.grad(2,2) + q.grad(3,1)))) + 6*q.val(4)*q.grad(3,2)*q.grad(4,1) + 4*sqrt(3)*q.val(4)*q.grad(1,1)*q.grad(4,1) + 6*q.val(5)*q.grad(3,2)*q.grad(5,1) + 4*sqrt(3)*q.val(5)*q.grad(1,1)*q.grad(5,1)));
    Fb = Fb + Ib;
    Fe = Fe + Ie;

    %Get Volume of Nematic Phase
    V = Real(1,1)
    V = V + Integral(Omega,(1/Sn.val)*(sqrt(3)*q.val(1) + 3*q.val(2)));
    
    Quadrature_Order = 4;
    
    G1 = GeoElement(Omega);
    
    MATS = Matrices(Quadrature_Order,G1);
    MATS = MATS.Append_Matrix(Fb);
    MATS = MATS.Append_Matrix(Fe);
    MATS = MATS.Append_Matrix(V);

end
