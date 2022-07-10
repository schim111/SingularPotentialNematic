function MATS = MatAssem2D_MS3()

    %define domain
    Omega = Domain('triangle');
    
    %finite element spaces
    Vector_P1_qn = Element(Omega,lagrange_deg1_dim2,5);
    Vector_P1_dq = Element(Omega,lagrange_deg1_dim2,3);
    
    %FE space functions
    vv = Test(Vector_P1_dq);
    uu = Trial(Vector_P1_dq);
    
    %Coef functions and constants
    q_n = Coef(Vector_P1_qn);
    H = Coef(Vector_P1_dq);
    alph = Constant(Omega);
    L = Constant(Omega);
    L2 = Constant(Omega);
    L3 = Constant(Omega);
    dt = Constant(Omega);
    
    %FE matrices
    Mass_Matrix1 = Bilinear(Vector_P1_dq,Vector_P1_dq);
    I1 = Integral(Omega,uu.val(1)*vv.val(1));
    I2 = Integral(Omega,uu.val(2)*vv.val(2));
    I3 = Integral(Omega,uu.val(3)*vv.val(3));
    Mass_Matrix1 = Mass_Matrix1 + I1 + I2 + I3;
    
    Stiff_Matrix1 = Bilinear(Vector_P1_dq,Vector_P1_dq);
    I4 = Integral(Omega,4*L.val*dt.val*(uu.grad(1)' * vv.grad(1)));
    I5 = Integral(Omega,4*L.val*dt.val*(uu.grad(2)' * vv.grad(2)));
    I6 = Integral(Omega,4*L.val*dt.val*(uu.grad(3)' * vv.grad(3)));
    I7 = Integral(Omega,2*L2.val*dt.val*(uu.grad(1)'*vv.grad(1)));
    I8 = Integral(Omega,2*L2.val*dt.val*uu.grad(2,1)*vv.grad(2,1));
    I9 = Integral(Omega,2*L2.val*dt.val*uu.grad(3,2)*vv.grad(2,1));
    I10 = Integral(Omega,2*L2.val*dt.val*uu.grad(2,1)*vv.grad(3,2));
    I11 = Integral(Omega,2*L2.val*dt.val*uu.grad(3,2)*vv.grad(3,2));
    I17 = Integral(Omega,L3.val*dt.val*(8/sqrt(3))*q_n.val(1)*uu.grad(1,1)*vv.grad(1,1));
    I18 = Integral(Omega,L3.val*dt.val*(4*q_n.val(2)*uu.grad(1,2) - (4/sqrt(3))*q_n.val(1)*uu.grad(1,2))*vv.grad(1,2));
    I19 = Integral(Omega,L3.val*dt.val*((8/sqrt(3))*q_n.val(1)*uu.grad(2,1))*vv.grad(2,1));
    I20 = Integral(Omega,L3.val*dt.val*(4*q_n.val(2)*uu.grad(2,2) - (4/sqrt(3))*q_n.val(1)*uu.grad(2,2))*vv.grad(2,2));
    I21 = Integral(Omega,L3.val*dt.val*((8/sqrt(3))*q_n.val(1)*uu.grad(3,1))*vv.grad(3,1));
    I22 = Integral(Omega,L3.val*dt.val*(4*q_n.val(2)*uu.grad(3,2) - (4/sqrt(3))*q_n.val(1)*uu.grad(3,2))*vv.grad(3,2));
    Stiff_Matrix1 = Stiff_Matrix1 + I4 + I5 + I6 + I7 + I8 + I9 + I10 + I11 + I17 + I18 + I19 + I20 + I21 + I22;

    RHS1 = Linear(Vector_P1_dq);
    I12 = Integral(Omega,(1 + 4*alph.val*dt.val)*q_n.val(3)*vv.val(1));
    I13 = Integral(Omega,(1 + 4*alph.val*dt.val)*q_n.val(4)*vv.val(2));
    I14 = Integral(Omega,(1 + 4*alph.val*dt.val)*q_n.val(5)*vv.val(3));
    I15 = Integral(Omega,-2*L2.val*dt.val*(-(1/sqrt(3))*q_n.grad(1,2) + q_n.grad(2,2))*vv.grad(1,1));
    I16 = Integral(Omega,-(4/sqrt(3))*L2.val*dt.val*q_n.grad(1,1)*vv.grad(1,2));
    I23 = Integral(Omega,L3.val*dt.val*(-4*q_n.grad(1,2)*q_n.grad(1,1) - 4*q_n.grad(2,2)*q_n.grad(2,1))*vv.val(1));
    I24 = Integral(Omega,2*dt.val*(H.val(1)*H.val(2)*vv.val(1) + H.val(1)*H.val(3)*vv.val(2) + H.val(2)*H.val(3)*vv.val(3)));
    RHS1 = RHS1 + I12 + I13 + I14 + I15 + I16 + I23 + I24;
    
    %Set minimum quad order
    Quadrature_Order = 4;
    %Define geometry representation
    G1 = GeoElement(Omega);
    %Define set of matrices
    MATS = Matrices(Quadrature_Order,G1);
    
    %collect MATS together
    MATS = MATS.Append_Matrix(Mass_Matrix1);
    MATS = MATS.Append_Matrix(Stiff_Matrix1);
    MATS = MATS.Append_Matrix(RHS1);
end
