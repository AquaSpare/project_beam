%%% Projection yttre clamped, olika inre proj/SAT mix

%%%% För exakt lösning IC=1
%%%% proj_projSAT2_proj(41,-1,0,1,0.4,1,100,4,1,500,0.0001, 1, 1, 1);

%%%% Några exempel:

%%%% Clamped:
%%%% proj_projSAT2_proj(81,-1,0,1,0.1,1,100,4,1,200, 0.00001, 0, 1, 1);

%%%% Free:
%%%% proj_projSAT2_proj(81,-1,0,1,0.1,1,100,4,2,200, 0.00001, 0, 1, 1);

%%%% Sliding
%%%% proj_projSAT2_proj(81,-1,0,1,0.1,1,100,4,3,200, 0.00001, 0, 1, 1);

%%%% Hinged
%%%% proj_projSAT2_proj(81,-1,0,1,0.1,1,100,4,4,200, 0.00001, 0, 1, 1);
function [w,t,x,h,k,error] = proj_projSAT2_proj(N, x0, xl, xN, T, a1, a2, order, BC, plotspeed, k_ratio, IC, plotornot,timestepperversion,timestepperOrder)
close all
pause on
%%%Domain%%%
h = (xN-xl)/(N-1);
x1 = x0:h:xl;
x2 = xl:h:xN;
x = [x1 x2];
%%%%%%%%%%%%%%%%%%%%%%%%

%%% Initial condition %%%
if(IC == 0)
    xr = -1/4;
    r0 = 1/30;  
    u = @(x) exp(-(xr-x).^2/r0^2);
    v = @(x) exp(-(xr-x).^2/r0^2);
    u0 = u(x1);
    v0 = v(x2);
else 
    [u,v] = beam_ana(BC, a1, a2);
    u0 = u(x1,0);
    v0 = v(x2,0);
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

w0 = [u0 v0];
w0_t = zeros(1,length(w0));

if timestepperversion == 1
    [w,k,t] = timestepper(T,h,A,w0,w0_t, k_ratio,timestepperOrder);
elseif timestepperversion == 2
    [w,k,t] = timestepperv2(T,h,A,w0,w0_t,k_ratio,timestepperOrder);
elseif timestepperversion == 3
    [w,k,t,error] = timestepperv3(T,h,A,w0,w0_t,k_ratio,BC,a1,a2,x0,xl,xN,timestepperOrder);
end

%%% Plot
if plotornot == 1
    figure(1);
    if(IC ~= 0)
        ymax = max(abs(u(x1,0)));
    end
    for i = 1:plotspeed:length(t)
        if(IC ~= 0)
            plot(x1, w(1:N,i), '*b', x2, w(N+1:end,i), '*r');
            hold on;
            plot(x1,u(x1,t(i)), 'b', x2, v(x2,t(i)), 'r');
            axis([x0 xN -ymax ymax]);
            hold off;
        else
            plot(x1, w(1:N,i), 'b', x2, w(N+1:end,i), 'r');
            axis([x0 xN -1 1]);
        end
        pause(0.00000001);
    end
end