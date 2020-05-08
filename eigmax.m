function maxeig = eigmax(method,N,ratio,BC,order,b1,b2)
%returns maximum eigenvalue of the iteration matrix for different methods

%METHOD == 1 for proj_proj_proj
%METHOD == 2 for proj_projSAT_proj
%METHOD == 3 for SAT_projSAT_SAT
%METHOD == 4 for homogenous beam with projection method
%METHOD == 5 for homogenous beam with SAT method

if method == 1 || method == 2 || method == 3
    h = 2/(N-1);
else 
    h = 1/(N-1);
end

k = ratio*h^2;

if (order == 2)
        [D4, H, HI, M, e1, eN, d1_1, dN_1, d1_2, dN_2, d1_3, dN_3] = SBP2_D4(N, h);
elseif (order == 4)
        [D4, H, HI, M, e1, eN, d1_1, dN_1, d1_2, dN_2, d1_3, dN_3] = SBP4_D4(N, h);
elseif (order == 6)
        [D4, H, HI, M, e1, eN, d1_1, dN_1, d1_2, dN_2, d1_3, dN_3] = SBP6_D4(N, h);
end

if method == 1    
    if(BC == 1) %clamped
        L = [eN -e1; dN_1 -d1_1; b1*dN_2 -b2*d1_2; b1*dN_3 -b2*d1_3; e1 zeros(1,N); d1_1 zeros(1,N); zeros(1,N) eN; zeros(1,N) dN_1];
    elseif(BC == 2) %free
        L = [eN -e1; dN_1 -d1_1; b1*dN_2 -b2*d1_2; b1*dN_3 -b2*d1_3; d1_2 zeros(1,N); d1_3 zeros(1,N); zeros(1,N) dN_2; zeros(1,N) dN_3];
    elseif(BC == 3) %sliding
        L = [eN -e1; dN_1 -d1_1; b1*dN_2 -b2*d1_2; b1*dN_3 -b2*d1_3; d1_1 zeros(1,N); d1_3 zeros(1,N); zeros(1,N) dN_1; zeros(1,N) dN_3];
    elseif(BC == 4) %hinged
        L = [eN -e1; dN_1 -d1_1; b1*dN_2 -b2*d1_2; b1*dN_3 -b2*d1_3; e1 zeros(1,N); d1_2 zeros(1,N); zeros(1,N) eN; zeros(1,N) dN_2];
    end
    
    H = [H zeros(N); zeros(N) H];
    HI = inv(H);
    P = [eye(N) zeros(N); zeros(N) eye(N)] - HI*L'*inv(L*HI*L')*L;
    A = (-P*[b1*D4 zeros(N); zeros(N) b2*D4]*P);
    
    maxeig = max(abs(eig(k^2.*A+2*eye(size(A)))));
end


if method == 2
    %%%Boundary Conditions%%%
    if(BC == 1)%clamped
        L = [eN -e1; dN_1 -d1_1; e1 zeros(1,N); d1_1 zeros(1,N); zeros(1,N) eN; zeros(1,N) dN_1];
    elseif(BC == 2)%free
        L = [eN -e1; dN_1 -d1_1; d1_2 zeros(1,N); d1_3 zeros(1,N); zeros(1,N) dN_2; zeros(1,N) dN_3];
    elseif(BC == 3)%sliding
        L = [eN -e1; dN_1 -d1_1; d1_1 zeros(1,N); d1_3 zeros(1,N); zeros(1,N) dN_1; zeros(1,N) dN_3];
    elseif(BC == 4)%hinged
        L = [eN -e1; dN_1 -d1_1; e1 zeros(1,N); d1_2 zeros(1,N); zeros(1,N) eN; zeros(1,N) dN_2];
    end
    l_u = -b1*D4 + b1*HI*(eN'*dN_3 - dN_1'*dN_2);
    r_l = -b2*D4 + b2*HI*(d1_1'*d1_2 - e1'*d1_3);

    H = [H zeros(N); zeros(N) H];
    HI = inv(H);
    P = [eye(N) zeros(N); zeros(N) eye(N)] - HI*L'*inv(L*HI*L')*L;

    A = P*[l_u zeros(N); zeros(N) r_l]*P;
    
    maxeig = max(abs(eig(k^2.*A+2*eye(size(A)))));
end

if method == 3
    
    %Bounary conditions
    if(BC == 2)%free
    L = [eN -e1; dN_1 -d1_1];
    l_u = -b1*D4 + b1*HI*(d1_1'*d1_2 - e1'*d1_3 + eN'*dN_3 - dN_1'*dN_2);
    r_l = -b2*D4 + b2*HI*(d1_1'*d1_2 - e1'*d1_3 + eN'*dN_3 - dN_1'*dN_2);
   elseif(BC == 1) %clamped
    tau3 = 1/2;
    tau4 = 1/2;
    tau5 = 1-tau3;
    tau6 = 1-tau4;
    L = [eN -e1; dN_1 -d1_1];
    if(order == 4)
        alfa_2 = 0.548;
        alfa_3 = 1.088;
    end
    tau1 = 4/(h^3*alfa_3);
    tau2 = 4/(h*alfa_2);
    tau7 = tau1;
    tau8 = tau2;
    
    l_u = b1*D4 +b1*HI*(-d1_3'*e1 + d1_2'*d1_1 + (tau1*e1)'*e1 +(tau2*d1_1)'*d1_1 -tau3*eN'*dN_3 + tau4*dN_1'*dN_2);
    r_u = b2*HI*(tau3*eN'*d1_3 - tau4*dN_1'*d1_2);

    l_l = b1*HI*(-tau5*e1'*dN_3 + tau6*d1_1'*dN_2);
    r_l = b2*D4 + b2*HI*(tau5*e1'*d1_3 - tau6*d1_1'*d1_2 +dN_3'*eN + (tau7*eN')*eN - dN_2'*dN_1 + (tau8*dN_1')*dN_1);

    end
    
    H = [H zeros(N); zeros(N) H];
    HI = inv(H);
    P = [eye(N) zeros(N); zeros(N) eye(N)] - HI*L'*inv(L*HI*L')*L;

    A = P*[l_u zeros(N); zeros(N) r_l]*P;
    maxeig= max(abs(eig(k^2.*A+2*eye(size(A)))));
end

if method == 4
    
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
    A = -b1.*P*D4*P;
    maxeig= max(abs(eig(k^2.*A+2*eye(size(A)))));
end



if method == 5 %SAT homogenous beam
    
   if order == 2
       alfa_2 = 1.250;
       alfa_3 = 0.4;
   elseif order == 4
       alfa_2 = 0.548;
       alfa_3 = 1.088;
   else
       alfa_2 = 0.322;
       alfa_3 = 0.156;
   end
   
   if BC == 1 %clamped
       tau0_1 = 4/(h^3*alfa_3);
       tauL_1 = tau0_1;
       tau0_2 = 4/(h*alfa_2);
       tauL_2 = tau0_2;
       SAT = HI*(d1_3-tau0_1*e1)'*e1 - HI*(d1_2+tau0_2*d1_1)'*d1_1 - HI*(dN_3+tauL_1*eN)'*eN + HI*(dN_2-tauL_2*dN_1)'*dN_1;
   elseif BC == 2 %free
       SAT = HI*d1_1'*d1_2 -HI*e1'*d1_3 -HI*dN_1'*dN_2 + HI*eN'*dN_3;
   end
   A = -D4 + SAT;
   maxeig= max(abs(eig(k^2.*A+2*eye(size(A)))));
end



end

