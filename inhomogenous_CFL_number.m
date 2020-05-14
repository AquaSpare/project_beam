function CFL_number = inhomogenous_CFL_number(method,N,BC,a1,a2)
%method 1, PURE PROJECTION
%method 2, 2nd, 3rd SAT
%method 3, 2nd SAT
%method 4, 3rd SAT
xN = 1;
xl = 0;
h = (xN-xl)/(N-1);
i = 1;
for order = 2:2:6
    
    if(order == 2)
        [D4, H, HI, M, e1, eN, d1_1, dN_1, d1_2, dN_2, d1_3, dN_3] = SBP2_D4(N, h);
    elseif(order == 4)
        [D4, H, HI, M, e1, eN, d1_1, dN_1, d1_2, dN_2, d1_3, dN_3] = SBP4_D4(N, h);
    elseif(order == 6)
        [D4, H, HI, M, e1, eN, d1_1, dN_1, d1_2, dN_2, d1_3, dN_3] = SBP6_D4(N, h);
    end
    
    if method == 1
                %%%Boundary Conditions%%%
        if(BC == 1) %clamped
            L = [eN -e1; dN_1 -d1_1; a1*dN_2 -a2*d1_2; a1*dN_3 -a2*d1_3; e1 zeros(1,N); d1_1 zeros(1,N); zeros(1,N) eN; zeros(1,N) dN_1];
        elseif(BC == 2) %free
            L = [eN -e1; dN_1 -d1_1; a1*dN_2 -a2*d1_2; a1*dN_3 -a2*d1_3; d1_2 zeros(1,N); d1_3 zeros(1,N); zeros(1,N) dN_2; zeros(1,N) dN_3];
        elseif(BC == 3) %sliding
            L = [eN -e1; dN_1 -d1_1; a1*dN_2 -a2*d1_2; a1*dN_3 -a2*d1_3; d1_1 zeros(1,N); d1_3 zeros(1,N); zeros(1,N) dN_1; zeros(1,N) dN_3];
        elseif(BC == 4) %hinged
            L = [eN -e1; dN_1 -d1_1; a1*dN_2 -a2*d1_2; a1*dN_3 -a2*d1_3; e1 zeros(1,N); d1_2 zeros(1,N); zeros(1,N) eN; zeros(1,N) dN_2];
        end

        H = [H zeros(N); zeros(N) H];
        HI = inv(H);
        P = [eye(N) zeros(N); zeros(N) eye(N)] - HI*L'*inv(L*HI*L')*L;
        A = (-P*[a1*D4 zeros(N); zeros(N) a2*D4]*P);

    elseif method == 2
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

        tau3 = 1/2;
        tau4 = 1/2;
        tau5 = 1-tau3;
        tau6 = 1-tau4;

        l_u = a1*D4 + a1*HI*(-tau3*eN'*dN_3 + tau4*dN_1'*dN_2);
        r_u = a2*HI*(tau3*eN'*d1_3 - tau4*dN_1'*d1_2);

        l_l = a1*HI*(-tau5*e1'*dN_3 + tau6*d1_1'*dN_2);
        r_l = a2*D4 + a2*HI*(tau5*e1'*d1_3 - tau6*d1_1'*d1_2);

        H = [H zeros(N); zeros(N) H];
        HI = inv(H);
        P = [eye(N) zeros(N); zeros(N) eye(N)] - HI*L'*inv(L*HI*L')*L;

        A = -P*[l_u r_u; l_l r_l]*P;
    elseif method == 3
                %%%Boundary Conditions%%%
        if(BC == 1)%clamped
            L = [eN -e1; dN_1 -d1_1; a1*dN_3 -a2*d1_3; e1 zeros(1,N); d1_1 zeros(1,N); zeros(1,N) eN; zeros(1,N) dN_1];
        elseif(BC == 2)%free
            L = [eN -e1; dN_1 -d1_1; a1*dN_3 -a2*d1_3; d1_2 zeros(1,N); d1_3 zeros(1,N); zeros(1,N) dN_2; zeros(1,N) dN_3];
        elseif(BC == 3)%sliding
            L = [eN -e1; dN_1 -d1_1; a1*dN_3 -a2*d1_3; d1_1 zeros(1,N); d1_3 zeros(1,N); zeros(1,N) dN_1; zeros(1,N) dN_3];
        elseif(BC == 4)%hinged
            L = [eN -e1; dN_1 -d1_1; a1*dN_3 -a2*d1_3; e1 zeros(1,N); d1_2 zeros(1,N); zeros(1,N) eN; zeros(1,N) dN_2];
        end

        tau4 = 1/2;
        tau6 = 1-tau4;

        l_u = a1*D4 + a1*HI*(tau4*dN_1'*dN_2);
        r_u = a2*HI*(- tau4*dN_1'*d1_2);

        l_l = a1*HI*(tau6*d1_1'*dN_2);
        r_l = a2*D4 + a2*HI*(- tau6*d1_1'*d1_2);

        H = [H zeros(N); zeros(N) H];
        HI = inv(H);
        P = [eye(N) zeros(N); zeros(N) eye(N)] - HI*L'*inv(L*HI*L')*L;

        A = -P*[l_u r_u; l_l r_l]*P;
        
    elseif method == 4
                %%%Boundary Conditions%%%
        if(BC == 1)%clamped
            L = [eN -e1; dN_1 -d1_1; a1*dN_2 -a2*d1_2; e1 zeros(1,N); d1_1 zeros(1,N); zeros(1,N) eN; zeros(1,N) dN_1];
        elseif(BC == 2)%free
            L = [eN -e1; dN_1 -d1_1; a1*dN_2 -a2*d1_2; d1_2 zeros(1,N); d1_3 zeros(1,N); zeros(1,N) dN_2; zeros(1,N) dN_3];
        elseif(BC == 3)%sliding
            L = [eN -e1; dN_1 -d1_1; a1*dN_2 -a2*d1_2; d1_1 zeros(1,N); d1_3 zeros(1,N); zeros(1,N) dN_1; zeros(1,N) dN_3];
        elseif(BC == 4)%hinged
            L = [eN -e1; dN_1 -d1_1; a1*dN_2 -a2*d1_2; e1 zeros(1,N); d1_2 zeros(1,N); zeros(1,N) eN; zeros(1,N) dN_2];
        end

        tau3 = 1/2;
        tau5 = 1-tau3;

        l_u = a1*D4 + a1*HI*(-tau3*eN'*dN_3);
        r_u = a2*HI*(tau3*eN'*d1_3);

        l_l = a1*HI*(-tau5*e1'*dN_3);
        r_l = a2*D4 + a2*HI*(tau5*e1'*d1_3);

        H = [H zeros(N); zeros(N) H];
        HI = inv(H);
        P = [eye(N) zeros(N); zeros(N) eye(N)] - HI*L'*inv(L*HI*L')*L;

        A = -P*[l_u r_u; l_l r_l]*P;
    end
    CFL_number(i) = 2/sqrt(max(abs(eig(h^4*A))));
    i=i+1;
end