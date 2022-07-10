function MATS = Get_Mu_Assemble()

    %Define Domain
    Omega = Domain('triangle');

    %Define finite element spaces
    Vector_P1 = Element(Omega,lagrange_deg1_dim2,5);
    
    %Functions
    q = Coef(Vector_P1);
    Lam = Coef(Vector_P1);
    alph = Constant(Omega);
    L = Constant(Omega);
    LS = Constant(Omega);
    L2 = Constant(Omega);
    L3 = Constant(Omega);
    L4 = Constant(Omega);
    

    %Define eval matrices
    mu = Real(1,1);
    I1 = Integral(Omega,sqrt(3)*(Lam.val(1) - 4*alph.val*q.val(1) - 4*L.val*(q.hess(1,1,1) + q.hess(2,2,1)) - LS.val*(6*(q.hess(1,1,1) + q.hess(2,2,1)) + 6*sqrt(3)*(q.hess(1,1,2) + q.hess(2,2,2))) - L2.val*((8/3)*q.hess(1,1,1) + (4/sqrt(3))*q.hess(1,2,3) + (2/3)*q.hess(2,2,1) - (2/sqrt(3))*q.hess(2,2,2) - (2/sqrt(3))*q.hess(2,1,3)) - L3.val*(2/sqrt(3))*(q.grad(1,2)^2 + q.grad(2,2)^2 + q.grad(3,2)^2 + q.grad(4,2)^2 + q.grad(5,2)^2 - 2*(q.grad(1,1)^2 + q.grad(2,1)^2 + q.grad(3,1)^2 + q.grad(4,1)^2 + q.grad(5,1)^2)) + L4.val*(-(2/sqrt(3))*q.grad(1,2)^2 + (4/sqrt(3))*q.grad(1,1)^2 + 2*q.grad(1,1)*q.grad(3,2) + 2*q.grad(1,2)*q.grad(2,2) + 2*q.grad(1,2)*q.grad(3,1))));
    I2 = Integral(Omega,3*(Lam.val(2) - 4*alph.val*q.val(2) - 4*L.val*(q.hess(1,1,2) + q.hess(2,2,2)) - LS.val*(18*(q.hess(1,1,2) + q.hess(2,2,2)) + 6*sqrt(3)*(q.hess(1,1,1) + q.hess(2,2,1))) - 2*L2.val*(-(1/sqrt(3))*q.hess(2,2,1) + q.hess(2,2,2) + q.hess(2,1,3)) + 2*L3.val*(q.grad(1,2)^2 + q.grad(2,2)^2 + q.grad(3,2)^2 + q.grad(4,2)^2 + q.grad(5,2)^2) + L4.val*(2*q.grad(2,1)*q.grad(3,2) + (4/sqrt(3))*q.grad(1,1)*q.grad(2,1) - (2/sqrt(3))*q.grad(1,2)*q.grad(2,2) + 2*q.grad(2,2)^2 + 2*q.grad(2,2)*q.grad(3,1))));
    mu = mu + I1 + I2;

    Quadrature_Order = 4;

    G1 = GeoElement(Omega);

    MATS = Matrices(Quadrature_Order,G1);
    MATS = MATS.Append_Matrix(mu);

end
