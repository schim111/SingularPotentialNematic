function MATS = MatAssem2D_MS4()

    %Define domain
    Omega = Domain('tetrahedron');
    
    %finite element spaces
    Vector_P1 = Element(Omega,lagrange_deg1_dim3,3);
    Vector_P1_dLam = Element(Omega,lagrange_deg1_dim3,9);
    Vector_P1_qn = Element(Omega,lagrange_deg1_dim3,5);
    
    %FE space functions
    vv = Test(Vector_P1);
    uu = Trial(Vector_P1);
    
    %Coef functions and constants
    V = Coef(Vector_P1)
    q_bar = Coef(Vector_P1);
    q_n = Coef(Vector_P1_qn);
    Lam = Coef(Vector_P1);
    dLam = Coef(Vector_P1_dLam);
    gamm = Constant(Omega);
    xi = Constant(Omega);
    L = Constant(Omega);
    L2 = Constant(Omega);
    L3 = Constant(Omega);
    dt = Constant(Omega);


    g = (1/sqrt(3))*q_n.val(1)*(2*V.grad(1,1) - V.grad(2,2) - V.grad(3,3)) + q_n.val(2)*(V.grad(2,2) - V.grad(3,3)) + q_bar.val(1)*(V.grad(2,1) + V.grad(1,2)) + q_bar.val(2)*(V.grad(3,1) + V.grad(1,3)) + q_bar.val(3)*(V.grad(3,2) + V.grad(2,3));
    
    %FE matrices 
    Mass_Matrix2 = Bilinear(Vector_P1,Vector_P1);
    I1 = Integral(Omega,gamm.val*dt.val*dLam.val(1)*uu.val(1)*vv.val(1));
    I2 = Integral(Omega,gamm.val*dt.val*dLam.val(2)*uu.val(2)*vv.val(1));
    I3 = Integral(Omega,gamm.val*dt.val*dLam.val(3)*uu.val(3)*vv.val(1));
    I4 = Integral(Omega,gamm.val*dt.val*dLam.val(4)*uu.val(1)*vv.val(2));
    I5 = Integral(Omega,gamm.val*dt.val*dLam.val(5)*uu.val(2)*vv.val(2));
    I6 = Integral(Omega,gamm.val*dt.val*dLam.val(6)*uu.val(3)*vv.val(2));
    I7 = Integral(Omega,gamm.val*dt.val*dLam.val(7)*uu.val(1)*vv.val(3));
    I8 = Integral(Omega,gamm.val*dt.val*dLam.val(8)*uu.val(2)*vv.val(3));
    I9 = Integral(Omega,gamm.val*dt.val*dLam.val(9)*uu.val(3)*vv.val(3));
    I19 = Integral(Omega,4*gamm.val*L3.val*dt.val*(q_bar.grad(1,2)*uu.grad(1,1) + q_bar.grad(1,1)*uu.grad(1,2) + q_bar.grad(2,2)*uu.grad(2,1) + q_bar.grad(2,1)*uu.grad(2,2) + q_bar.grad(3,2)*uu.grad(3,1) + q_bar.grad(3,1)*uu.grad(3,2))*vv.val(1));
    I20 = Integral(Omega,4*gamm.val*L3.val*dt.val*(q_bar.grad(1,1)*uu.grad(1,3) + q_bar.grad(1,3)*uu.grad(1,1) + q_bar.grad(2,1)*uu.grad(2,3) + q_bar.grad(2,3)*uu.grad(2,1) + q_bar.grad(3,1)*uu.grad(3,3) + q_bar.grad(3,3)*uu.grad(3,1))*vv.val(2));
    I21 = Integral(Omega,4*gamm.val*L3.val*dt.val*(q_bar.grad(1,2)*uu.grad(1,3) + q_bar.grad(1,3)*uu.grad(1,2) + q_bar.grad(2,2)*uu.grad(2,3) + q_bar.grad(2,3)*uu.grad(2,2) + q_bar.grad(3,2)*uu.grad(3,3) + q_bar.grad(3,3)*uu.grad(3,2))*vv.val(3));
    IV4 = Integral(Omega,dt.val*xi.val*(2*uu.val(1)*g + 2*q_bar.val(1)*uu.val(1)*(V.grad(2,1) + V.grad(1,2)))*vv.val(1));
    IV5 = Integral(Omega,dt.val*xi.val*(2*uu.val(2)*g + 2*q_bar.val(2)*uu.val(2)*(V.grad(3,1) + V.grad(1,3)))*vv.val(2));
    IV6 = Integral(Omega,dt.val*xi.val*(2*uu.val(3)*g + 2*q_bar.val(3)*uu.val(3)*(V.grad(3,2) + V.grad(2,3)))*vv.val(3));
    Mass_Matrix2 = Mass_Matrix2 + I1 + I2 + I3 + I4 + I5 + I6 + I7 + I8 + I9 + I19 + I20 + I21 + IV4 + IV5 + IV6;

    Stiff_Matrix2 = Bilinear(Vector_P1,Vector_P1);
    I22 = Integral(Omega,4*gamm.val*L3.val*dt.val*(q_bar.grad(1,2)*uu.val(1) + q_bar.val(1)*uu.grad(1,2) + q_bar.grad(1,3)*uu.val(2) + q_bar.val(2)*uu.grad(1,3))*vv.grad(1,1));
    I23 = Integral(Omega,4*gamm.val*L3.val*dt.val*(q_bar.grad(1,1)*uu.val(1) + q_bar.val(1)*uu.grad(1,1) + q_bar.grad(1,3)*uu.val(3) + q_bar.val(3)*uu.grad(1,3))*vv.grad(1,2));
    I24 = Integral(Omega,4*gamm.val*L3.val*dt.val*(q_bar.grad(1,1)*uu.val(2) + q_bar.val(2)*uu.grad(1,1) + q_bar.grad(1,2)*uu.val(3) + q_bar.val(3)*uu.grad(1,2))*vv.grad(1,3));
    I25 = Integral(Omega,4*gamm.val*L3.val*dt.val*(q_bar.grad(2,2)*uu.val(1) + q_bar.val(1)*uu.grad(2,2) + q_bar.grad(2,3)*uu.val(2) + q_bar.val(2)*uu.grad(2,3))*vv.grad(2,1));
    I26 = Integral(Omega,4*gamm.val*L3.val*dt.val*(q_bar.grad(2,1)*uu.val(1) + q_bar.val(1)*uu.grad(2,1) + q_bar.grad(2,3)*uu.val(3) + q_bar.val(3)*uu.grad(2,3))*vv.grad(2,2));
    I27 = Integral(Omega,4*gamm.val*L3.val*dt.val*(q_bar.grad(2,1)*uu.val(2) + q_bar.val(2)*uu.grad(2,1) + q_bar.grad(2,2)*uu.val(3) + q_bar.val(3)*uu.grad(2,2))*vv.grad(2,3));
    I28 = Integral(Omega,4*gamm.val*L3.val*dt.val*(q_bar.grad(3,2)*uu.val(1) + q_bar.val(1)*uu.grad(3,2) + q_bar.grad(3,3)*uu.val(2) + q_bar.val(2)*uu.grad(3,3))*vv.grad(3,1));
    I29 = Integral(Omega,4*gamm.val*L3.val*dt.val*(q_bar.grad(3,1)*uu.val(1) + q_bar.val(1)*uu.grad(3,1) + q_bar.grad(3,3)*uu.val(3) + q_bar.val(3)*uu.grad(3,3))*vv.grad(3,2));
    I30 = Integral(Omega,4*gamm.val*L3.val*dt.val*(q_bar.grad(3,1)*uu.val(2) + q_bar.val(2)*uu.grad(3,1) + q_bar.grad(3,2)*uu.val(3) + q_bar.val(3)*uu.grad(3,2))*vv.grad(3,3));
    Stiff_Matrix2 = Stiff_Matrix2  + I22 + I23 + I24 + I25 + I26 + I27 + I28 + I29 + I30;
    
    RHS2 = Linear(Vector_P1);
    I10 = Integral(Omega,(-q_bar.val(1) - gamm.val*dt.val*Lam.val(1))*vv.val(1));
    I11 = Integral(Omega,-4*gamm.val*L.val*dt.val*(q_bar.grad(1)' * vv.grad(1)));
    I12 = Integral(Omega,(-q_bar.val(2) - gamm.val*dt.val*Lam.val(2))*vv.val(2));
    I13 = Integral(Omega,-4*gamm.val*L.val*dt.val*(q_bar.grad(2)' * vv.grad(2)));
    I14 = Integral(Omega,(-q_bar.val(3) - gamm.val*dt.val*Lam.val(3))*vv.val(3));
    I15 = Integral(Omega,-4*gamm.val*L.val*dt.val*(q_bar.grad(3)' * vv.grad(3)));
    I16 = Integral(Omega,-2*gamm.val*L2.val*dt.val*(q_bar.grad(1,1)*vv.grad(1,1) + q_bar.grad(3,3)*vv.grad(1,1) + q_bar.grad(1,2)*vv.grad(1,2) + q_bar.grad(2,3)*vv.grad(1,2)));
    I17 = Integral(Omega,-2*gamm.val*L2.val*dt.val*(q_bar.grad(2,1)*vv.grad(2,1) + q_bar.grad(3,2)*vv.grad(2,1) + q_bar.grad(1,2)*vv.grad(2,3) + q_bar.grad(2,3)*vv.grad(2,3)));
    I18 = Integral(Omega,-2*gamm.val*L2.val*dt.val*(q_bar.grad(2,1)*vv.grad(3,2) + q_bar.grad(3,2)*vv.grad(3,2) + q_bar.grad(1,1)*vv.grad(3,3) + q_bar.grad(3,3)*vv.grad(3,3)));
    I31 = Integral(Omega,-4*gamm.val*L3.val*dt.val*(q_bar.grad(1,2)*q_bar.grad(1,1) + q_bar.grad(2,2)*q_bar.grad(2,1) + q_bar.grad(3,2)*q_bar.grad(3,1))*vv.val(1));
    I32 = Integral(Omega,gamm.val*L3.val*dt.val*(-4*q_bar.val(1)*q_bar.grad(1,2) - (8/sqrt(3))*q_n.val(1)*q_bar.grad(1,1) - 4*q_bar.val(2)*q_bar.grad(1,3))*vv.grad(1,1));
    I33 = Integral(Omega,gamm.val*L3.val*dt.val*(-4*q_n.val(2)*q_bar.grad(1,2) - 4*q_bar.val(1)*q_bar.grad(1,1) + (4/sqrt(3))*q_n.val(1)*q_bar.grad(1,2) - 4*q_bar.val(3)*q_bar.grad(1,3))*vv.grad(1,2));
    I34 = Integral(Omega,gamm.val*L3.val*dt.val*((4/sqrt(3))*q_n.val(1)*q_bar.grad(1,3) + 4*q_n.val(2)*q_bar.grad(1,3) - 4*q_bar.val(2)*q_bar.grad(1,1) - 4*q_bar.val(3)*q_bar.grad(1,2))*vv.grad(1,3));
    I35 = Integral(Omega,-4*gamm.val*L3.val*dt.val*(q_bar.grad(1,1)*q_bar.grad(1,3) + q_bar.grad(2,1)*q_bar.grad(2,3) + q_bar.grad(3,1)*q_bar.grad(3,3))*vv.val(2));
    I36 = Integral(Omega,gamm.val*L3.val*dt.val*(-4*q_bar.val(1)*q_bar.grad(2,2) - (8/sqrt(3))*q_n.val(1)*q_bar.grad(2,1) - 4*q_bar.val(2)*q_bar.grad(2,3))*vv.grad(2,1));
    I37 = Integral(Omega,gamm.val*L3.val*dt.val*(-4*q_n.val(2)*q_bar.grad(2,2) - 4*q_bar.val(1)*q_bar.grad(2,1) + (4/sqrt(3))*q_n.val(1)*q_bar.grad(2,2) - 4*q_bar.val(3)*q_bar.grad(2,3))*vv.grad(2,2));
    I38 = Integral(Omega,gamm.val*L3.val*dt.val*((4/sqrt(3))*q_n.val(1)*q_bar.grad(2,3) + 4*q_n.val(2)*q_bar.grad(2,3) - 4*q_bar.val(2)*q_bar.grad(2,1) - 4*q_bar.val(3)*q_bar.grad(2,2))*vv.grad(2,3));
    I39 = Integral(Omega,-4*gamm.val*L3.val*dt.val*(q_bar.grad(1,2)*q_bar.grad(1,3) + q_bar.grad(2,2)*q_bar.grad(2,3) + q_bar.grad(3,2)*q_bar.grad(3,3))*vv.val(3));
    I40 = Integral(Omega,gamm.val*L3.val*dt.val*(-4*q_bar.val(1)*q_bar.grad(3,2) - (8/sqrt(3))*q_n.val(1)*q_bar.grad(3,1) - 4*q_bar.val(2)*q_bar.grad(3,3))*vv.grad(3,1));
    I41 = Integral(Omega,gamm.val*L3.val*dt.val*(-4*q_n.val(2)*q_bar.grad(3,2) - 4*q_bar.val(1)*q_bar.grad(3,1) + (4/sqrt(3))*q_n.val(1)*q_bar.grad(3,2) - 4*q_bar.val(3)*q_bar.grad(3,3))*vv.grad(3,2));
    I42 = Integral(Omega,gamm.val*L3.val*dt.val*((4/sqrt(3))*q_n.val(1)*q_bar.grad(3,3) + 4*q_n.val(2)*q_bar.grad(3,3) - 4*q_bar.val(2)*q_bar.grad(3,1) - 4*q_bar.val(3)*q_bar.grad(3,2))*vv.grad(3,3));
    IV1 = Integral(Omega,-dt.val*((V.val(1)*q_bar.grad(1,1) + V.val(2)*q_bar.grad(1,2) + V.val(3)*q_bar.grad(1,3)) + xi.val*(2*q_bar.val(1)*g - q_bar.val(1)*(V.grad(1,1) + V.grad(2,2)) - 0.5*q_bar.val(2)*(V.grad(3,2) + V.grad(2,3)) - 0.5*q_bar.val(3)*(V.grad(3,1) + V.grad(1,3))) - 0.5*(q_bar.val(2)*(V.grad(3,2) - V.grad(2,3)) - q_bar.val(3)*(V.grad(3,1) - V.grad(1,3))))*vv.val(1));
    IV2 = Integral(Omega,-dt.val*((V.val(1)*q_bar.grad(2,1) + V.val(2)*q_bar.grad(2,2) + V.val(3)*q_bar.grad(2,3)) + xi.val*(2*q_bar.val(2)*g - 0.5*q_bar.val(1)*(V.grad(3,2) + V.grad(2,3)) - q_bar.val(2)*(V.grad(1,1) + V.grad(3,3)) - 0.5*q_bar.val(3)*(V.grad(2,1) + V.grad(1,2))) - 0.5*(q_bar.val(1)*(V.grad(3,2) - V.grad(2,3)) - q_bar.val(3)*(V.grad(2,1) - V.grad(1,2))))*vv.val(2));
    IV3 = Integral(Omega,-dt.val*((V.val(1)*q_bar.grad(3,1) + V.val(2)*q_bar.grad(3,2) + V.val(3)*q_bar.grad(3,3)) + xi.val*(2*q_bar.val(3)*g - 0.5*q_bar.val(1)*(V.grad(3,1) + V.grad(1,3)) - 0.5*q_bar.val(2)*(V.grad(2,1) + V.grad(1,2) - q_bar.val(3)*(V.grad(2,2) + V.grad(3,3)))) - 0.5*(q_bar.val(1)*(V.grad(3,1) - V.grad(1,3)) - q_bar.val(2)*(V.grad(1,2) - V.grad(2,1))))*vv.val(3));
    RHS2 = RHS2 + I10 + I11 + I12 + I13 + I14 + I15 + I16 + I17 + I18 + I31 + I32 + I33 + I34 + I35 + I36 + I37 + I38 + I39 + I40 + I41 + I42 + IV1 + IV2 + IV3;
    
    %Set minimum quad order
    Quadrature_Order = 4;
    %Define Geometry representation
    G1 = GeoElement(Omega);
    %Define set of matrices
    MATS = Matrices(Quadrature_Order,G1);
    
    %collect Mats together
    MATS = MATS.Append_Matrix(Mass_Matrix2);
    MATS = MATS.Append_Matrix(RHS2);
    MATS = MATS.Append_Matrix(Stiff_Matrix2);
    
end
