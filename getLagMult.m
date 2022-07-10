%Newton's method for Lagrange multipliers
%Returns multipliers for every input q in same format (q should be Mx5)
%Also returns gradient of Lamda wrt q
function [Lam, dLam1, dLam2] = getLagMult(q,tol,LEB)
    Lam = 0*q;
    M = length(q(:,1));
    dLam1 = zeros(M,4);
    dLam2 = zeros(M,9);
    
    
    %Go through each deg. of freedom
    parfor ii=1:M
        for newt=1:30 %Newton iterations
            %Get Z RHS and Jacobian
            Z = ZFunc(Lam(ii,:),LEB);
            GradZ = ZGrad(Lam(ii,:),LEB);
            HessZ = ZHess(Lam(ii,:),LEB);
            RHS = -(1/Z)*GradZ + [q(ii,1) + (sqrt(3)/6);q(ii,2) + 0.5;q(ii,3);q(ii,4);q(ii,5)];
            Jac = (1/Z)*HessZ - (1/Z^2)*(GradZ * GradZ');
            
            %update
            deltaLam = Jac \ RHS;
            
            Lam(ii,:) = Lam(ii,:) + deltaLam';
            
            %Check convergence
            err = max(abs(deltaLam));
            if err < tol
                break;
            end
            if newt==30
                disp('Newton Iteration did not converge!');
            end
        end
        
        %find dLam1 and dLam2
        dLam = (Jac \ eye(5));
        dLam = dLam'; %rows as columns
        dLamtemp1 = dLam(1:2,1:2);
        dLamtemp2 = dLam(3:5,3:5);
        dLam1(ii,:) = dLamtemp1(:)';
        dLam2(ii,:) = dLamtemp2(:)';
    end
end
