function [w,t,x,h,k] = beam_eq_projection_dav(N,x0,xl,xN,t0,T,b1,b2,order,BC,k_ratio,IC)
close all
pause on
%%%Domain%%%
h = (xN-xl)/(N-1);
x1 = x0:h:xl;
x2 = xl:h:xN;
x = [x1 x2];

%%%%%%%%%%%%%

%%%Constant%%%

%%%%%%%%%%%

%%%Initial data%%%
u0 = IC{1}(x1,0);
v0 = IC{2}(x2,0);

%%%SBP operators%%%
if(order == 2)
    [D4, H, HI, M, e1, eN, d1_1, dN_1, d1_2, dN_2, d1_3, dN_3] = SBP2_D4(N, h);
elseif(order == 4)
    [D4, H, HI, M, e1, eN, d1_1, dN_1, d1_2, dN_2, d1_3, dN_3] = SBP4_D4(N, h);
elseif(order == 6)
    [D4, H, HI, M, e1, eN, d1_1, dN_1, d1_2, dN_2, d1_3, dN_3] = SBP6_D4(N, h);
end

%%%Boundary Conditions%%%
if(BC == 1)%clamped
    L = [eN -e1; dN_1 -d1_1; b1*dN_2 -b2*d1_2; b1*dN_3 -b2*d1_3; e1 zeros(1,N); d1_1 zeros(1,N); zeros(1,N) eN; zeros(1,N) dN_1];
elseif(BC == 2)%free
    L = [eN -e1; dN_1 -d1_1; b1*dN_2 -b2*d1_2; b1*dN_3 -b2*d1_3; d1_2 zeros(1,N); d1_3 zeros(1,N); zeros(1,N) dN_2; zeros(1,N) dN_3];
elseif(BC == 3)%sliding
    L = [eN -e1; dN_1 -d1_1; b1*dN_2 -b2*d1_2; b1*dN_3 -b2*d1_3; d1_1 zeros(1,N); d1_3 zeros(1,N); zeros(1,N) dN_1; zeros(1,N) dN_3];
elseif(BC == 4)%hinged
    L = [eN -e1; dN_1 -d1_1; b1*dN_2 -b2*d1_2; b1*dN_3 -b2*d1_3; e1 zeros(1,N); d1_2 zeros(1,N); zeros(1,N) eN; zeros(1,N) dN_2];
end
H = [H zeros(N); zeros(N) H];
HI = inv(H);

P = [eye(N) zeros(N); zeros(N) eye(N)] - HI*L'*inv(L*HI*L')*L;

A = (-P*[b1^2*D4 zeros(N); zeros(N) b2^2*D4]*P);

w0 = [u0 v0];
w0_t = 0;

[w,k,t] = timestepper(T,h,A,w0,w0_t,k_ratio);
%%%Exact Solution%%%

