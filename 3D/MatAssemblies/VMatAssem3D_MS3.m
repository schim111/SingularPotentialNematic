function MATS = MatAssem2D_MS3()

    %define domain
    Omega = Domain('tetrahedron');
    
    %finite element spaces
    Vector_P1_qn = Element(Omega,lagrange_deg1_dim3,5);
    Vector_P1_dq = Element(Omega,lagrange_deg1_dim3,3);
    
    %FE space functions
    vv = Test(Vector_P1_dq);
    uu = Trial(Vector_P1_dq);
    
    %Coef functions and constants
    V = Coef(Vector_P1_dq);
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
    I3 = Integral(Omega,uu.val(3)*vv.val(3));
    IV1 = Integral(Omega,dt.val*(V.val(1)*uu.grad(1,1) + V.val(2)*uu.grad(1,2) + V.val(3)*uu.grad(1,3))*vv.val(1));
    IV2 = Integral(Omega,dt.val*(0.5*(V.grad(2,3) - V.grad(3,2))*uu.val(2)*vv.val(1)));
    IV3 = Integral(Omega,dt.val*0.5*(V.grad(3,1) - V.grad(1,3))*uu.val(3)*vv.val(1));
    IV4 = Integral(Omega,dt.val*(-0.5*(V.grad(3,2) - V.grad(2,3))*uu.val(1)*vv.val(2)));
    IV5 = Integral(Omega,dt.val*(V.val(1)*uu.grad(2,1) + V.val(2)*uu.grad(2,2) + V.val(3)*uu.grad(2,3))*vv.val(2));
    IV6 = Integral(Omega,dt.val*0.5*(V.grad(2,1) - V.grad(1,2))*uu.val(3)*vv.val(2));
    IV7 = Integral(Omega,dt.val*(-0.5*(V.grad(3,1) - V.grad(1,3))*uu.val(1)*vv.val(3)));
    IV8 = Integral(Omega,dt.val*0.5*(V.grad(1,2) - V.grad(2,1))*uu.val(2)*vv.val(3));
    IV9 = Integral(Omega,dt.val*(V.val(1)*uu.grad(3,1) + V.val(2)*uu.grad(3,2) + V.val(3)*uu.grad(3,3))*vv.val(3));
    IV13 = Integral(Omega,dt.val*xi.val*uu.val(1)*(V.grad(1,1) + V.grad(2,2))*vv.val(1));
    IV14 = Integral(Omega,-dt.val*xi.val*0.5*uu.val(2)*(V.grad(3,2) + V.grad(2,3))*vv.val(1));
    IV15 = Integral(Omega,-dt.val*xi.val*0.5*uu.val(3)*(V.grad(3,1) + V.grad(1,3))*vv.val(1));
    IV16 = Integral(Omega,-dt.val*xi.val*0.5*uu.val(1)*(V.grad(3,2) + V.grad(2,3))*vv.val(2));
    IV17 = Integral(Omega,-dt.val*xi.val*uu.val(2)*(V.grad(1,1) + V.grad(3,3))*vv.val(2));
    IV18 = Integral(Omega,-dt.val*xi.val*uu.val(3)*(V.grad(2,1) + V.grad(1,2))*vv.val(2));
    IV19 = Integral(Omega,-dt.val*xi.val*uu.val(1)*(V.grad(3,1) + V.grad(1,3))*vv.val(3));
    IV20 = Integral(Omega,-dt.val*xi.val*uu.val(2)*(V.grad(2,1) + V.grad(1,2))*vv.val(3));
    IV21 = Integral(Omega,-dt.val*xi.val*uu.val(3)*(V.grad(2,2) + V.grad(3,3))*vv.val(3)); 
    Mass_Matrix1 = Mass_Matrix1 + I1 + I2 + I3 + IV1 + IV2 + IV3 + IV4 + IV5 + IV6 + IV7 + IV8 + IV9 + IV13 + IV14 + IV15 + IV16 + IV17 + IV18 + IV19 + IV20 + IV21;
    
    Stiff_Matrix1 = Bilinear(Vector_P1_dq,Vector_P1_dq);
    I4 = Integral(Omega,4*gamm.val*L.val*dt.val*(uu.grad(1)' * vv.grad(1)));
    I5 = Integral(Omega,4*gamm.val*L.val*dt.val*(uu.grad(2)' * vv.grad(2)));
    I6 = Integral(Omega,4*gamm.val*L.val*dt.val*(uu.grad(3)' * vv.grad(3)));
    I7 = Integral(Omega,2*gamm.val*L2.val*dt.val*(uu.grad(1,1)*vv.grad(1,1) + uu.grad(3,3)*vv.grad(1,1) + uu.grad(1,2)*vv.grad(1,2) + uu.grad(2,3)*vv.grad(1,2)));
    I8 = Integral(Omega,2*gamm.val*L2.val*dt.val*(uu.grad(2,1)*vv.grad(2,1) + uu.grad(3,2)*vv.grad(2,1) + uu.grad(1,2)*vv.grad(2,3) + uu.grad(2,3)*vv.grad(2,3)));
    I9 = Integral(Omega,2*gamm.val*L2.val*dt.val*(uu.grad(2,1)*vv.grad(3,2) + uu.grad(3,2)*vv.grad(3,2) + uu.grad(1,1)*vv.grad(3,3) + uu.grad(3,3)*vv.grad(3,3)));
    I19 = Integral(Omega,gamm.val*L3.val*dt.val*(8/sqrt(3))*q_n.val(1)*uu.grad(1,1)*vv.grad(1,1));
    I20 = Integral(Omega,gamm.val*L3.val*dt.val*(4*q_n.val(2)*uu.grad(1,2) - (4/sqrt(3))*q_n.val(1)*uu.grad(1,2))*vv.grad(1,2));
    I21 = Integral(Omega,gamm.val*L3.val*dt.val*(-(4/sqrt(3))*q_n.val(1)*uu.grad(1,3) - 4*q_n.val(2)*uu.grad(1,3))*vv.grad(1,3));
    I22 = Integral(Omega,gamm.val*L3.val*dt.val*((8/sqrt(3))*q_n.val(1)*uu.grad(2,1))*vv.grad(2,1));
    I23 = Integral(Omega,gamm.val*L3.val*dt.val*(4*q_n.val(2)*uu.grad(2,2) - (4/sqrt(3))*q_n.val(1)*uu.grad(2,2))*vv.grad(2,2));
    I24 = Integral(Omega,gamm.val*L3.val*dt.val*(-(4/sqrt(3))*q_n.val(1)*uu.grad(2,3) - 4*q_n.val(2)*uu.grad(2,3))*vv.grad(2,3));
    I25 = Integral(Omega,gamm.val*L3.val*dt.val*((8/sqrt(3))*q_n.val(1)*uu.grad(3,1))*vv.grad(3,1));
    I26 = Integral(Omega,gamm.val*L3.val*dt.val*(4*q_n.val(2)*uu.grad(3,2) - (4/sqrt(3))*q_n.val(1)*uu.grad(3,2))*vv.grad(3,2));
    I27 = Integral(Omega,gamm.val*L3.val*dt.val*(-(4/sqrt(3))*q_n.val(1)*uu.grad(3,3) - 4*q_n.val(2)*uu.grad(3,3))*vv.grad(3,3));
    Stiff_Matrix1 = Stiff_Matrix1 + I4 + I5 + I6 + I7 + I8 + I9 + I19 + I20 + I21 + I22 + I23 + I24 + I25 + I26 + I27;

    RHS1 = Linear(Vector_P1_dq);
    I10 = Integral(Omega,(1 + 4*gamm.val*alph.val*dt.val)*q_n.val(3)*vv.val(1));
    I11 = Integral(Omega,(1 + 4*gamm.val*alph.val*dt.val)*q_n.val(4)*vv.val(2));
    I12 = Integral(Omega,(1 + 4*gamm.val*alph.val*dt.val)*q_n.val(5)*vv.val(3));
    I13 = Integral(Omega,-2*gamm.val*L2.val*dt.val*(-(1/sqrt(3))*q_n.grad(1,2) + q_n.grad(2,2))*vv.grad(1,1));
    I14 = Integral(Omega,-(4/sqrt(3))*gamm.val*L2.val*dt.val*q_n.grad(1,1)*vv.grad(1,2));
    I15 = Integral(Omega,-2*gamm.val*L2.val*dt.val*(-(1/sqrt(3))*q_n.grad(1,3) - q_n.grad(2,3))*vv.grad(2,1));
    I16 = Integral(Omega,-2*gamm.val*L2.val*dt.val*(2/sqrt(3))*q_n.grad(1,1)*vv.grad(2,3));
    I17 = Integral(Omega,-2*gamm.val*L2.val*dt.val*(-(1/sqrt(3))*q_n.grad(1,3) - q_n.grad(2,3))*vv.grad(3,2));
    I18 = Integral(Omega,-2*gamm.val*L2.val*dt.val*(-(1/sqrt(3))*q_n.grad(1,2) + q_n.grad(2,2))*vv.grad(3,3));
    I28 = Integral(Omega,-4*gamm.val*L3.val*dt.val*(q_n.grad(1,2)*q_n.grad(1,1) + q_n.grad(2,2)*q_n.grad(2,1))*vv.val(1));
    I29 = Integral(Omega,-4*gamm.val*L3.val*dt.val*(q_n.grad(1,1)*q_n.grad(1,3) + q_n.grad(2,1)*q_n.grad(2,3))*vv.val(2));
    I30 = Integral(Omega,-4*gamm.val*L3.val*dt.val*(q_n.grad(1,2)*q_n.grad(1,3) + q_n.grad(2,2)*q_n.grad(2,3))*vv.val(3));
    IV10 = Integral(Omega,-dt.val*(xi.val*(-0.5*(V.grad(2,1) + V.grad(1,2))*((1/sqrt(3))*q_n.val(1) + q_n.val(2)) - (1/3)*(V.grad(2,1) + V.grad(1,2))) - 0.5*(sqrt(3)*q_n.val(1) - q_n.val(2))*(V.grad(2,1) - V.grad(1,2)))*vv.val(1));
    IV11 = Integral(Omega,-dt.val*(xi.val*(-0.5*(V.grad(3,1) + V.grad(1,3))*((1/sqrt(3))*q_n.val(1) - q_n.val(2)) - (1/3)*(V.grad(3,1) + V.grad(1,3))) - 0.5*(sqrt(3)*q_n.val(1) + q_n.val(2))*(V.grad(3,1) - V.grad(1,3)))*vv.val(2));
    IV12 = Integral(Omega,-dt.val*(xi.val*((1/sqrt(3))*q_n.val(1)*(V.grad(3,2) + V.grad(2,3)) - (1/3)*(V.grad(3,2) + V.grad(2,3))) - q_n.val(2)*(V.grad(3,2) - V.grad(2,3)))*vv.val(3));
    RHS1 = RHS1 + I10 + I11 + I12 + I13 + I14 + I15 + I16 + I17 + I18 + I28 + I29 + I30 + IV10 + IV11 + IV12;
    
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
