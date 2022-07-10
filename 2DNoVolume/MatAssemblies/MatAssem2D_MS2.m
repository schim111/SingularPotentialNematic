function MATS = MatAssem2D_MS2()

    %define domain
    Omega = Domain('triangle');
    
    %finite element spaces
    Vector_P1 = Element(Omega,lagrange_deg1_dim2,2);
    Vector_P1_dLam = Element(Omega,lagrange_deg1_dim2,4);
    Vector_P1_qn = Element(Omega,lagrange_deg1_dim2,5);
    
    %FE space functions
    vv = Test(Vector_P1);
    uu = Trial(Vector_P1);

    %Coef functions and constants
    q_bar = Coef(Vector_P1);
    q_n = Coef(Vector_P1_qn);
    Lam = Coef(Vector_P1);
    dLam = Coef(Vector_P1_dLam);
    L = Constant(Omega);
    L2 = Constant(Omega);
    L3 = Constant(Omega);
    dt = Constant(Omega);
    
    %FE matrices
    Mass_Matrix2 = Bilinear(Vector_P1,Vector_P1);
    I1 = Integral(Omega,dt.val*dLam.val(1)*uu.val(1)*vv.val(1));
    I2 = Integral(Omega,dt.val*dLam.val(2)*uu.val(2)*vv.val(1));
    I3 = Integral(Omega,dt.val*dLam.val(3)*uu.val(1)*vv.val(2));
    I4 = Integral(Omega,dt.val*dLam.val(4)*uu.val(2)*vv.val(2));
    I12 = Integral(Omega,L3.val*dt.val*(-(4/sqrt(3))*q_bar.grad(1,2)*uu.grad(1,2) + (8/sqrt(3))*q_bar.grad(1,1)*uu.grad(1,1) - (4/sqrt(3))*q_bar.grad(2,2)*uu.grad(2,2) + (8/sqrt(3))*q_bar.grad(2,1)*uu.grad(2,1))*vv.val(1));
    I13 = Integral(Omega,L3.val*dt.val*(4*q_bar.grad(1,2)*uu.grad(1,2) + 4*q_bar.grad(2,2)*uu.grad(2,2))*vv.val(2));
    Mass_Matrix2 = Mass_Matrix2 + I1 + I2 + I3 + I4 + I12 + I13;

    Stiff_Matrix2 = Bilinear(Vector_P1,Vector_P1);
    I14 = Integral(Omega,L3.val*dt.val*((8/sqrt(3))*q_bar.grad(1,1)*uu.val(1) + (8/sqrt(3))*q_bar.val(1)*uu.grad(1,1))*vv.grad(1,1));
    I15 = Integral(Omega,L3.val*dt.val*(4*q_bar.val(2)*uu.grad(1,2) - (4/sqrt(3))*q_bar.grad(1,2)*uu.val(1) - (4/sqrt(3))*q_bar.val(1)*uu.grad(1,2) + 4*q_bar.grad(1,2)*uu.val(2))*vv.grad(1,2));
    I16 = Integral(Omega,L3.val*dt.val*((8/sqrt(3))*q_bar.grad(2,1)*uu.val(1) + (8/sqrt(3))*q_bar.val(1)*uu.grad(2,1))*vv.grad(2,1));
    I17 = Integral(Omega,L3.val*dt.val*(-(4/sqrt(3))*q_bar.grad(2,2)*uu.val(1) + 4*q_bar.grad(2,2)*uu.val(2) + 4*q_bar.val(2)*uu.grad(2,2) - (4/sqrt(3))*q_bar.val(1)*uu.grad(2,2))*vv.grad(2,2));
    Stiff_Matrix2 = Stiff_Matrix2 + I14 + I15 + I16 + I17;
    
    RHS2 = Linear(Vector_P1);
    I5 = Integral(Omega,(-q_bar.val(1) - dt.val*Lam.val(1))*vv.val(1));
    I6 = Integral(Omega,-4*L.val*dt.val*(q_bar.grad(1)' * vv.grad(1)));
    I7 = Integral(Omega,(-q_bar.val(2) - dt.val*Lam.val(2))*vv.val(2));
    I8 = Integral(Omega,-4*L.val*dt.val*(q_bar.grad(2)' * vv.grad(2)));
    I9 = Integral(Omega,-(8/3)*L2.val*dt.val*q_bar.grad(1,1)*vv.grad(1,1));
    I10 = Integral(Omega,2*L2.val*dt.val*(-(1/3)*q_bar.grad(1,2) + (1/sqrt(3))*q_bar.grad(2,2))*vv.grad(1,2));
    I11 = Integral(Omega,-2*L2.val*dt.val*(-(1/sqrt(3))*q_bar.grad(1,2) + q_bar.grad(2,2))*vv.grad(2,2));
    I18 = Integral(Omega,L3.val*dt.val*(1/sqrt(3))*(2*q_bar.grad(1,2)^2 + 2*q_bar.grad(2,2)^2 - 4*q_bar.grad(1,1)^2 - 4*q_bar.grad(2,1)^2)*vv.val(1));
    I19 = Integral(Omega,L3.val*dt.val*(-4*q_n.val(3)*q_bar.grad(1,2) - (8/sqrt(3))*q_bar.val(1)*q_bar.grad(1,1))*vv.grad(1,1));
    I20 = Integral(Omega,L3.val*dt.val*(-4*q_bar.val(2)*q_bar.grad(1,2) - 4*q_n.val(3)*q_bar.grad(1,1) + (4/sqrt(3))*q_bar.val(1)*q_bar.grad(1,2))*vv.grad(1,2));
    I21 = Integral(Omega,L3.val*dt.val*(-2*q_bar.grad(1,2)^2 - 2*q_bar.grad(2,2)^2)*vv.val(2));
    I22 = Integral(Omega,L3.val*dt.val*(-4*q_n.val(3)*q_bar.grad(2,2) - (8/sqrt(3))*q_bar.val(1)*q_bar.grad(2,1))*vv.grad(2,1));
    I23 = Integral(Omega,L3.val*dt.val*(-4*q_bar.val(2)*q_bar.grad(2,2) - 4*q_n.val(3)*q_bar.grad(2,1) + (4/sqrt(3))*q_bar.val(1)*q_bar.grad(2,2))*vv.grad(2,2));
    RHS2 = RHS2 + I5 + I6 + I7 + I8 + I9 + I10 + I11 + I18 + I19 + I20 + I21 + I22 + I23;
    
    %Set minimum quad order
    Quadrature_Order = 4;
    %Define geometry representation
    G1 = GeoElement(Omega);
    %Define set of matrices
    MATS = Matrices(Quadrature_Order,G1);
    
    %Collects Mats together
    MATS = MATS.Append_Matrix(Mass_Matrix2);
    MATS = MATS.Append_Matrix(RHS2);
    MATS = MATS.Append_Matrix(Stiff_Matrix2);
    
end
