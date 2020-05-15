function [w, t, x, h, k, error] = SAT3_enBalk(N,x0,xN,T, a, order, BC, plotspeed, k_ratio, IC, plotornot, timestepperversion)

h = (xN-x0)/(N-1);
x = x0:h:xN;

%%% Initial condition %%%
if(IC == 0)
    xr = 1/4;
    r0 = 1/30;
    u = @(x) exp(-(xr-x).^2/r0^2);
    u0 = u(x);
    u0_t = 0;
else 
    u = homo_beam_ana_2(a, BC);
    u0 = u(x,0);
    u0_t = 0;
end

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
    SAT = HI*((d1_3-tau0_1*e1)'*e1 - (dN_3+tauL_1*eN)'*eN);
    L = [d1_1;dN_1];
elseif BC == 2 %free
    SAT = HI*(d1_1'*d1_2 -dN_1'*dN_2);
    L = [d1_3;dN_3];
elseif BC == 3 %sliding
    tau = 2/(h*alfa_2);
    SAT = HI*(-(d1_2+ tau*d1_1)'*d1_1 + (dN_2-tau*dN_1)'*dN_1);
    L = [d1_3;dN_3];
elseif BC == 4 %hinged
    tau = 4/(h^3 * alfa_3);
    SAT = HI*((d1_3-tau*e1)'*e1 -(dN_3+tau*eN)'*eN);
    L = [d1_2;dN_2];
end

P = eye(N)-HI*L'*inv(L*HI*L')*L;
A = a.*P*(-D4 + SAT)*P;
if timestepperversion == 1
    [w,k,t] = timestepper(T,h,A,u0,u0_t,k_ratio);
    error = 0;
elseif timestepperversion == 2
    [w,k,t] = timestepperv2(T,h,A,u0,u0_t,k_ratio);
    error = 0;
elseif timestepperversion == 3
    [w,k,t,error] = timestepperv3_enBalk(T,h,A,u0,u0_t,k_ratio,BC,a,x0,xN);
end



if plotornot == 1
    figure(1);
    for i = 1:plotspeed:length(t)
        if(IC == 1)
          plot(x, w(1:N,i), '*b');
           hold on;
           plot(x,u(x,t(i)), 'g');
           axis([0 1 -2 2]);
           hold off;
        else 
         plot(x, w(1:N,i));
         axis([x0 xN -1 1]);
        end
     pause(0.00000001);
    end
    
end


