function [u, t, x, h, k] = SAT_enBalk(N,x0,xN,T, b, order, BC, plotspeed, k_ratio,plotornot, timestepperversion)

h = (xN-x0)/(N-1);
x = x0:h:xN;

%%%Exact Solution and initial conditions
[u_exact,u0,u0_t] = homo_beam_ana(x,BC);

%%%%%%%

%%%SBP operators%%%
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

HI = inv(H);

if BC == 1 %clamped
    tau0_1 = 4/(h^3*alfa_3);
    tauL_1 = tau0_1;
    tau0_2 = 4/(h*alfa_2);
    tauL_2 = tau0_2;
    SAT = HI*(d1_3-tau0_1*e1)'*e1 - HI*(d1_2+tau0_2*d1_1)'*d1_1 - HI*(dN_3+tauL_1*eN)'*eN + HI*(dN_2-tauL_2*dN_1)'*dN_1;
elseif BC == 2 %free
    SAT = HI*d1_1'*d1_2 -HI*e1'*d1_3 -HI*dN_1'*dN_2 + HI*eN'*dN_3;
end

A = b.*(-D4 + SAT);
if timestepperversion == 1
    [u,k,t] = timestepper(T,h,A,u0,u0_t,k_ratio);
elseif timestepperversion == 2
    [u,k,t] = timestepperv2(T,h,A,u0,u0_t,k_ratio);
else
    return
end


if plotornot == 1
    figure(1);
    for i = 1:plotspeed:length(t)
        if(x0 == 0 && xN == 1)
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
end


