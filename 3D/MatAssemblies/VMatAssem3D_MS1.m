function MATS = MatAssem2D_MS1()

    %define domain
    Omega = Domain('tetrahedron');
    
    %finite element spaces
    Vector_P1_qn = Element(Omega,lagrange_deg1_dim3,5);
    Vector_P1_dq = Element(Omega,lagrange_deg1_dim3,2);
    Vector_P1_V = Element(Omega,lagrange_deg1_dim3,3);
    
    %FE space functions
    vv = Test(Vector_P1_dq);
    uu = Trial(Vector_P1_dq);
   
    %Coef functions and constants
    V = Coef(Vector_P1_V);
    q_n = Coef(Vector_P1_qn);
    alph = Constant(Omega);
    gamm = Constant(Omega);
    xi = Constant(Omega);
    L = Constant(Omega);
    L2 = Constant(Omega);
    L3 = Constant(Omega);
    dt = Constant(Omega);
    
    %FE matrices
    Mass_Matrix1 = Bilinear(Vector_P1_dq,Vector_P1_dq);
    I1 = Integral(Omega,uu.val(1)*vv.val(1));
    I2 = Integral(Omega,uu.val(2)*vv.val(2));
    IV1 = Integral(Omega,dt.val*(V.val(1)*uu.grad(1,1) + V.val(2)*uu.grad(1,2) + V.val(3)*uu.grad(1,3))*vv.val(1));
    IV2 = Integral(Omega,dt.val*(V.val(1)*uu.grad(2,1) + V.val(2)*uu.grad(2,2) + V.val(3)*uu.grad(2,3))*vv.val(2));
    IV5 = Integral(Omega,dt.val*(-2*V.grad(1,1)*uu.val(1))*vv.val(1));
    IV6 = Integral(Omega,dt.val*((1/sqrt(3))*(V.grad(2,2)*uu.val(1) - V.grad(3,3)*uu.val(1)))*vv.val(2));
    IV7 = Integral(Omega,dt.val*(-V.grad(2,2)*uu.val(2) - V.grad(2,2)*uu.val(2))*vv.val(2)); 
    Mass_Matrix1 = Mass_Matrix1 + I1 + I2 + IV1 + IV2 + IV6 + IV7;
    
    Stiff_Matrix1 = Bilinear(Vector_P1_dq,Vector_P1_dq);
    I3 = Integral(Omega,(4*gamm.val*L.val*dt.val*(uu.grad(1)' * vv.grad(1)) + 4*L.val*dt.val*(uu.grad(2)' * vv.grad(2))));
    I4 = Integral(Omega,(8/3)*gamm.val*L2.val*dt.val*uu.grad(1,1)*vv.grad(1,1));
    I5 = Integral(Omega,(2/3)*gamm.val*L2.val*dt.val*uu.grad(1,2)*vv.grad(1,2));
    I6 = Integral(Omega,(2/3)*gamm.val*L2.val*dt.val*uu.grad(1,3)*vv.grad(1,3));
    I7 = Integral(Omega,-(2/sqrt(3))*gamm.val*L2.val*dt.val*uu.grad(2,2)*vv.grad(1,2));
    I8 = Integral(Omega,(2/sqrt(3))*gamm.val*L2.val*dt.val*uu.grad(2,3)*vv.grad(1,3));
    I9 = Integral(Omega,-(2/sqrt(3))*gamm.val*L2.val*dt.val*uu.grad(1,2)*vv.grad(2,2));
    I10 = Integral(Omega,2*gamm.val*L2.val*dt.val*uu.grad(2,2)*vv.grad(2,2));
    I11 = Integral(Omega,(2/sqrt(3))*gamm.val*L2.val*dt.val*uu.grad(1,3)*vv.grad(2,3));
    I12 = Integral(Omega,2*gamm.val*L2.val*dt.val*uu.grad(2,3)*vv.grad(2,3));
    I25 = Integral(Omega,4*gamm.val*L3.val*dt.val*(q_n.val(3)*uu.grad(1,2) + q_n.val(4)*uu.grad(1,3))*vv.grad(1,1));
    I26 = Integral(Omega,4*gamm.val*L3.val*dt.val*(q_n.val(3)*uu.grad(1,1) + q_n.val(5)*uu.grad(1,3))*vv.grad(1,2));
    I27 = Integral(Omega,4*gamm.val*L3.val*dt.val*(q_n.val(4)*uu.grad(1,1) + q_n.val(5)*uu.grad(1,2))*vv.grad(1,3));
    I28 = Integral(Omega,4*gamm.val*L3.val*dt.val*(q_n.val(3)*uu.grad(2,2) + q_n.val(4)*uu.grad(2,3))*vv.grad(2,1));
    I29 = Integral(Omega,4*gamm.val*L3.val*dt.val*(q_n.val(3)*uu.grad(2,1) + q_n.val(5)*uu.grad(2,3))*vv.grad(2,2));
    I30 = Integral(Omega,4*gamm.val*L3.val*dt.val*(q_n.val(4)*uu.grad(2,1) + q_n.val(5)*uu.grad(2,2))*vv.grad(2,3));
    Stiff_Matrix1 = Stiff_Matrix1 + I3 + I4 + I5 + I6 + I7 + I8 + I9 + I10 + I11 + I12 + I25 + I26 + I27 + I28 + I29 + I30;
    
    RHS1 = Linear(Vector_P1_dq);
    I13 = Integral(Omega,(1 + 4*gamm.val*alph.val*dt.val)*q_n.val(1)*vv.val(1));
    I14 = Integral(Omega,(1 + 4*gamm.val*alph.val*dt.val)*q_n.val(2)*vv.val(2));
    I15 = Integral(Omega,-(4/sqrt(3))*gamm.val*L2.val*dt.val*q_n.grad(3,2)*vv.grad(1,1));
    I16 = Integral(Omega,-(4/sqrt(3))*gamm.val*L2.val*dt.val*q_n.grad(4,3)*vv.grad(1,1));
    I17 = Integral(Omega,(2/sqrt(3))*gamm.val*L2.val*dt.val*q_n.grad(3,1)*vv.grad(1,2));
    I18 = Integral(Omega,(2/sqrt(3))*gamm.val*L2.val*dt.val*q_n.grad(5,3)*vv.grad(1,2));
    I19 = Integral(Omega,(2/sqrt(3))*gamm.val*L2.val*dt.val*q_n.grad(4,1)*vv.grad(1,3));
    I20 = Integral(Omega,(2/sqrt(3))*gamm.val*L2.val*dt.val*q_n.grad(5,2)*vv.grad(1,3));
    I21 = Integral(Omega,-2*gamm.val*L2.val*dt.val*q_n.grad(3,1)*vv.grad(2,2));
    I22 = Integral(Omega,-2*gamm.val*L2.val*dt.val*q_n.grad(5,3)*vv.grad(2,2));
    I23 = Integral(Omega,2*gamm.val*L2.val*dt.val*q_n.grad(4,1)*vv.grad(2,3));
    I24 = Integral(Omega,2*gamm.val*L2.val*dt.val*q_n.grad(5,2)*vv.grad(2,3));
    I31 = Integral(Omega,gamm.val*L3.val*dt.val*(1/sqrt(3))*(2*(q_n.grad(3,2)^2 + q_n.grad(4,2)^2 + q_n.grad(5,2)^2 + q_n.grad(3,3)^2 + q_n.grad(4,3)^2 + q_n.grad(5,3)^2) - 4*(q_n.grad(3,1)^2 + q_n.grad(4,1)^2 + q_n.grad(5,1)^2))*vv.val(1));
    I32 = Integral(Omega,gamm.val*L3.val*dt.val*(-2*(q_n.grad(3,2)^2 + q_n.grad(4,2)^2 + q_n.grad(5,2)^2 - q_n.grad(3,3)^2 - q_n.grad(4,3)^2 - q_n.grad(5,3)^2))*vv.val(2));
    IV3 = Integral(Omega,-dt.val*(xi.val*(-0.5*sqrt(3)*(V.grad(2,1) + V.grad(1,2))*q_n.val(3) + 0.5*sqrt(3)*(V.grad(3,1) + V.grad(1,3))*q_n.val(4) + (1/sqrt(3))*V.grad(1,1)) - 0.5*sqrt(3)*(q_n.val(3)*(V.grad(1,2) - V.grad(2,1)) + q_n.val(4)*(V.grad(1,3) - V.grad(3,1))))*vv.val(1));
    IV4 = Integral(Omega,-dt.val*(xi.val*(-0.5*(V.grad(2,1) + V.grad(1,2))*q_n.val(3) + 0.5*(V.grad(3,1) + V.grad(1,3))*q_n.val(4) - (1/3)*(V.grad(2,2) - V.grad(3,3))) - 0.5*(q_n.val(3)*(V.grad(2,1) - V.grad(1,2)) - q_n.val(4)*(V.grad(3,1) - V.grad(1,3)) + 2*q_n.val(5)*(V.grad(2,3) - V.grad(3,2))))*vv.val(2));
    RHS1 = RHS1 + I13 + I14 + I15 + I16 + I17 + I18 + I19 + I20 + I21 + I22 + I23 + I24 + I31 + I32 + IV3 + IV4;
    
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
