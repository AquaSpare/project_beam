function [u, t, x, h, k] = SAT_enBalk(N,x0,xN,T, b, order, BC, plotspeed, k_ratio)

h = (xN-x0)/(N-1);
x = x0:h:xN;
%%%initial conditions%%%
%Fria ränder
u0 = @(x) cosh(1.50562*pi.*x) + cos(1.50562*pi.*x) - ((cos(1.50562*pi) -cosh(1.50562*pi))/(sin(1.50562*pi) - sinh(1.50562*pi)))*(sin(1.50562*pi.*x) +sinh(1.50562*pi.*x));
u0_t = 0;
%%%%%%%

%%%SBP operators%%%
if(order == 2)
    [D4, H, HI, M, e1, eN, d1_1, dN_1, d1_2, dN_2, d1_3, dN_3] = SBP2_D4(N, h);
elseif(order == 4)
    [D4, H, HI, M, e1, eN, d1_1, dN_1, d1_2, dN_2, d1_3, dN_3] = SBP4_D4(N, h);
elseif(order == 6)
    [D4, H, HI, M, e1, eN, d1_1, dN_1, d1_2, dN_2, d1_3, dN_3] = SBP6_D4(N, h);
end

if BC == 2 %free
    SAT = HI*d1_1'*d1_2 -HI*e1'*d1_3 -HI*dN_1'*dN_2 + HI*eN'*dN_3;
end

A = -D4 + SAT;

[u,k,t] = timestepper(T,h,A,u0(x),u0_t,k_ratio);
%%%Exact Solution%%%
u_exact = @(x,t) real(exp(-1i*(22.3733)*t)*u0(x));
figure(1);
for i = 1:plotspeed:length(t)
    if(BC == 2 && x0 == 0 && xN == 1)
        plot(x, u(1:N,i), '*b');
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

