function maxeig = bigboy(method,N,ratio,BC,order,b1,b2)

%METHOD == 1 for proj_proj_proj
%METHOD == 2 for proj_projSAT_proj
%METHOD == 3 for SAT_projSAT_SAT


h = 2/(N-1);
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
    else
        return;
    end
    
    H = [H zeros(N); zeros(N) H];
    HI = inv(H);
    P = [eye(N) zeros(N); zeros(N) eye(N)] - HI*L'*inv(L*HI*L')*L;

    A = P*[l_u zeros(N); zeros(N) r_l]*P;
    maxeig= max(abs(eig(k^2.*A+2*eye(size(A)))));
end


end

