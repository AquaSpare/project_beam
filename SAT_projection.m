
%%%% BARA BC = 2 %%%%

function [w,t,x,h,k] = SAT_projection(N, x0, xl, xN, T, a1, a2, order, BC, plotspeed, k_ratio, IC, plotornot, timestepperversion)
close all
pause on
%%%Domain%%%
h = (xN-xl)/(N-1);
x1 = x0:h:xl;
x2 = xl:h:xN;
x = [x1 x2];
%%%%%%%%%%%%%%%%%%%%%%%%

%%% Initial condition %%%
if(a1 == a2 && BC == 2 && x0 == 0 && xN == 1) %%%%%%exact solution exists 
    u0 = @(x) cosh(1.50562*pi.*x) + cos(1.50562*pi.*x) - ((cos(1.50562*pi) -cosh(1.50562*pi))/(sin(1.50562*pi) - sinh(1.50562*pi)))*(sin(1.50562*pi.*x) +sinh(1.50562*pi.*x));
    u0_t = 0;
else 
    u0 = IC{1}(x1,0);
    v0 = IC{2}(x2,0);
    u0_t = 0;
end

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
    L = [eN -e1; dN_1 -d1_1; e1 zeros(1,N); d1_1 zeros(1,N); zeros(1,N) eN; zeros(1,N) dN_1];
    %r_l = -a2*D4 + a2*HI*((d1_3 - e1'*d1_3 + eN'*dN_3 - dN_1'*dN_2); 
    %l_u = -a1*D4 + a1*HI*(d1_1'*d1_2 - e1'*d1_3 + eN'*dN_3 - dN_1'*dN_2);
elseif(BC == 2)%free
    L = [eN -e1; dN_1 -d1_1];
    l_u = -a1*D4 + a1*HI*(d1_1'*d1_2 - e1'*d1_3 + eN'*dN_3 - dN_1'*dN_2);
    r_l = -a2*D4 + a2*HI*(d1_1'*d1_2 - e1'*d1_3 + eN'*dN_3 - dN_1'*dN_2);
elseif(BC == 3)%sliding
    L = [eN -e1; dN_1 -d1_1; d1_1 zeros(1,N); d1_3 zeros(1,N); zeros(1,N) dN_1; zeros(1,N) dN_3];
    %r_l = -a2*D4 + a2*HI*(d1_1'*d1_2 - e1'*d1_3 + eN'*dN_3 - dN_1'*dN_2);
    %l_u = -a1*D4 + a1*HI*(d1_1'*d1_2 - e1'*d1_3 + eN'*dN_3 - dN_1'*dN_2);
elseif(BC == 4)%hinged
    L = [eN -e1; dN_1 -d1_1; e1 zeros(1,N); d1_2 zeros(1,N); zeros(1,N) eN; zeros(1,N) dN_2];
    %r_l = -a2*D4 + a2*HI*(d1_1'*d1_2 - e1'*d1_3 + eN'*dN_3 - dN_1'*dN_2);
    %l_u = -a1*D4 + a1*HI*(d1_1'*d1_2 - e1'*d1_3 + eN'*dN_3 - dN_1'*dN_2);
end

H = [H zeros(N); zeros(N) H];
HI = inv(H);
P = [eye(N) zeros(N); zeros(N) eye(N)] - HI*L'*inv(L*HI*L')*L;

A = P*[l_u zeros(N); zeros(N) r_l]*P;


w0 = [u0 v0];
w0_t = u0_t;

if timestepperversion == 1
    [w,k,t] = timestepper(T,h,A,w0,w0_t, k_ratio);
elseif timestepperversion == 2
    [w,k,t] = timestepperv2(T,h,A,w0,w0_t,k_ratio);
end

%%%Exact Solution%%%
u_exact = @(x,t) real(exp(-1i*(22.3733)*t)*u0(x));
if plotornot == 1
    figure(1);
    for i = 1:plotspeed:length(t)
        if(a1 == a2 && BC == 2 && x0 == 0 && xN == 1)
            plot(x1, w(1:N,i), '*b', x2, w(N+1:end,i), '*r');
            hold on;
            plot(x,u_exact(x,t(i)), 'g');
            axis([0 1 -2 2]);
            hold off;
        else 
            plot(x1, w(1:N,i), 'b', x2, w(N+1:end,i), 'r');
            axis([x0 xN -1 1]);
        end
        pause(0.00000001);
    end
end