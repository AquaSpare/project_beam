function CFL_number = homogenous_CFL_number(method,N,BC,a)
%method 1, PROJECTION
%method 2, SAT
h = 1/(N-1);
i=1;
for order = 2:2:6
    
if(order == 2)
    [D4, H, HI, M, e1, eN, d1_1, dN_1, d1_2, dN_2, d1_3, dN_3] = SBP2_D4(N, h);
    alfa_2 = 1.250;
    alfa_3 = 0.4;
elseif(order == 4)
    [D4, H, HI, M, e1, eN, d1_1, dN_1, d1_2, dN_2, d1_3, dN_3] = SBP4_D4(N, h);
    alfa_2 = 0.548;
    alfa_3 = 1.088;
elseif(order == 6)
    [D4, H, HI, M, e1, eN, d1_1, dN_1, d1_2, dN_2, d1_3, dN_3] = SBP6_D4(N, h);
    alfa_2 = 0.322;
    alfa_3 = 0.156;
end

if method == 1
        %------------Boundary conditions--------------
    if BC == 1
        %------Clamped-------
        L = [e1;eN;d1_1;dN_1];
    elseif BC == 2
        %-------Free---------
        L = [d1_2;dN_2;d1_3;dN_3];
    elseif BC == 3
        %------Sliding-------
        L = [d1_1;dN_1;d1_3;dN_3];
    elseif BC == 4
        %------Hinged--------
        L = [e1;eN;d1_2;dN_2];
    end
    
    
P = eye(N)-HI*L'*inv(L*HI*L')*L;
A = -a.*P*D4*P;

CFL_number(i) = 2/sqrt(max(abs(eig(h^4*A))));

elseif method == 2
    if BC == 1 %clamped
    tau0_1 = 4/(h^3*alfa_3);
    tauL_1 = tau0_1;
    tau0_2 = 4/(h*alfa_2);
    tauL_2 = tau0_2;
    SAT = HI*(d1_3-tau0_1*e1)'*e1 - HI*(d1_2+tau0_2*d1_1)'*d1_1 - HI*(dN_3+tauL_1*eN)'*eN + HI*(dN_2-tauL_2*dN_1)'*dN_1;
elseif BC == 2 %free
    SAT = HI*d1_1'*d1_2 -HI*e1'*d1_3 -HI*dN_1'*dN_2 + HI*eN'*dN_3;
elseif BC == 3 %sliding
    tau = 2/(h*alfa_2);
    SAT = HI*(-(d1_2+ tau*d1_1)'*d1_1 - e1'*d1_3 + (dN_2-tau*dN_1)'*dN_1 + eN'*dN_3);
elseif BC == 4 %hinged
    tau = 4/(h^3 * alfa_3);
    SAT = HI*((d1_3-tau*e1)'*e1 +d1_1'*d1_2 -(dN_3+tau*eN)'*eN -dN_1'*dN_2);
end

A = a.*(-D4 + SAT);

CFL_number(i) = 2/sqrt(max(abs(eig(h^4*A))));
end
i = i+1;
end