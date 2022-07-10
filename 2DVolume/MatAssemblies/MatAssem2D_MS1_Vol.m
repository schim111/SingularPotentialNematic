function MATS = MatAssem2D_MS1()

    %define domain
    Omega = Domain('triangle');
    
    %finite element spaces
    Vector_P1_qn = Element(Omega,lagrange_deg1_dim2,5);
    Vector_P1_dq = Element(Omega,lagrange_deg1_dim2,2);
    
    %FE space functions
    vv = Test(Vector_P1_dq);
    uu = Trial(Vector_P1_dq);
   
    %Coef functions and constants
    q_n = Coef(Vector_P1_qn);
    alph = Constant(Omega);
    L = Constant(Omega);
    LS = Constant(Omega);
    L2 = Constant(Omega);
    L3 = Constant(Omega);
    L4 = Constant(Omega);
    dt = Constant(Omega);
    mu = Constant(Omega);
    Sn = Constant(Omega);
    
    %FE matrices
    Mass_Matrix1 = Bilinear(Vector_P1_dq,Vector_P1_dq);
    I1 = Integral(Omega,uu.val(1)*vv.val(1));
    I2 = Integral(Omega,uu.val(2)*vv.val(2));
    I20 = Integral(Omega,L4.val*dt.val*(2*q_n.grad(3,2)*uu.grad(1,1) + 2*q_n.grad(3,1)*uu.grad(1,2))*vv.val(1));
    I21 = Integral(Omega,L4.val*dt.val*(2*q_n.grad(3,2)*uu.grad(2,1) + 2*q_n.grad(3,1)*uu.grad(2,2))*vv.val(2));
    Mass_Matrix1 = Mass_Matrix1 + I1 + I2 + I20 + I21;
    
    Stiff_Matrix1 = Bilinear(Vector_P1_dq,Vector_P1_dq);
    I3 = Integral(Omega,(4*L.val*dt.val*(uu.grad(1)' * vv.grad(1)) + 4*L.val*dt.val*(uu.grad(2)' * vv.grad(2))));
    I26 = Integral(Omega,LS.val*dt.val*((6*uu.grad(1)' + 6*sqrt(3)*uu.grad(2)') * vv.grad(1) + (18*uu.grad(2)' + 6*sqrt(3)*uu.grad(1)') * vv.grad(2)));
    I4 = Integral(Omega,(8/3)*L2.val*dt.val*uu.grad(1,1)*vv.grad(1,1));
    I5 = Integral(Omega,(2/3)*L2.val*dt.val*uu.grad(1,2)*vv.grad(1,2));
    I6 = Integral(Omega,-(2/sqrt(3))*L2.val*dt.val*uu.grad(2,2)*vv.grad(1,2));
    I7 = Integral(Omega,-(2/sqrt(3))*L2.val*dt.val*uu.grad(1,2)*vv.grad(2,2));
    I8 = Integral(Omega,2*L2.val*dt.val*uu.grad(2,2)*vv.grad(2,2));
    I14 = Integral(Omega,4*L3.val*dt.val*q_n.val(3)*uu.grad(1,2)*vv.grad(1,1));
    I15 = Integral(Omega,4*L3.val*dt.val*q_n.val(3)*uu.grad(1,1)*vv.grad(1,2));
    I16 = Integral(Omega,4*L3.val*dt.val*q_n.val(3)*uu.grad(2,2)*vv.grad(2,1));
    I17 = Integral(Omega,4*L3.val*dt.val*q_n.val(3)*uu.grad(2,1)*vv.grad(2,2));
    I22 = Integral(Omega,L4.val*dt.val*(2*q_n.grad(3,2)*uu.val(1))*vv.grad(1,1));
    I23 = Integral(Omega,L4.val*dt.val*(2*q_n.grad(3,1)*uu.val(2)*vv.grad(2,2) + 2*q_n.grad(3,2)*uu.val(2)*vv.grad(2,1)));
    Stiff_Matrix1 = Stiff_Matrix1 + I3 + I4 + I5 + I6 + I7 + I8 + I14 + I15 + I16 + I17 + I22 + I23 + I26;
    
    RHS1 = Linear(Vector_P1_dq);
    I9 = Integral(Omega,((1 + 4*alph.val*dt.val)*q_n.val(1) + dt.val*(sqrt(3)*mu.val/Sn.val))*vv.val(1));
    I10 = Integral(Omega,((1 + 4*alph.val*dt.val)*q_n.val(2) + dt.val*(3*mu.val/Sn.val))*vv.val(2));
    I11 = Integral(Omega,-(4/sqrt(3))*L2.val*dt.val*q_n.grad(3,2)*vv.grad(1,1));
    I12 = Integral(Omega,(2/sqrt(3))*L2.val*dt.val*q_n.grad(3,1)*vv.grad(1,2));
    I13 = Integral(Omega,-2*L2.val*dt.val*q_n.grad(3,1)*vv.grad(2,2));
    I18 = Integral(Omega,L3.val*dt.val*(1/sqrt(3))*(2*(q_n.grad(3,2)^2 + q_n.grad(4,2)^2 + q_n.grad(5,2)^2) - 4*(q_n.grad(3,1)^2 + q_n.grad(4,1)^2 + q_n.grad(5,1)^2))*vv.val(1));
    I19 = Integral(Omega,L3.val*dt.val*(-2*(q_n.grad(3,2)^2 + q_n.grad(4,2)^2 + q_n.grad(5,2)^2))*vv.val(2));
    I24 = Integral(Omega,L4.val*dt.val*((-4*sqrt(3)/3)*(q_n.val(3)*q_n.grad(3,1) + q_n.val(4)*q_n.grad(4,1) + q_n.val(5)*q_n.grad(5,1))*vv.grad(1,1) + (2*sqrt(3)/3)*(q_n.val(3)*q_n.grad(3,2) + q_n.val(4)*q_n.grad(4,2) + q_n.val(5)*q_n.grad(5,2))*vv.grad(1,2)));
    I25 = Integral(Omega,L4.val*dt.val*(-2*(q_n.val(3)*q_n.grad(3,2) + q_n.val(4)*q_n.grad(4,2) + q_n.val(5)*q_n.grad(5,2))*vv.grad(2,2)));
    RHS1 = RHS1 + I9 + I10 + I11 + I12 + I13 + I18 + I19 + I24 + I25;
    
    %Set minimum quad order
    Quadrature_Order = 4;
    %Define geometry representation
    G1 = GeoElement(Omega);
    %Define set of matrices
    MATS = Matrices(Quadrature_Order,G1);
    
    %Collect Mats together
    MATS = MATS.Append_Matrix(Mass_Matrix1);
    MATS = MATS.Append_Matrix(Stiff_Matrix1);
    MATS = MATS.Append_Matrix(RHS1);

end
