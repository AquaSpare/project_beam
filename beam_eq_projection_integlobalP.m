%%%% För exakt lösning krävs x0=0 xN=1, b1=b2=1 och bc=2 i.e
%%%% beam_eq_projection(81,0,1/2,1,0.5,1,1,4,2,100,0.001);

%%%% Några exempel:

%%%% Clamped:
%%%% beam_eq_projection(81,-1,0,1,0.1,1,100,4,1,200, 0.00001);

%%%% Free:
%%%% beam_eq_projection(81,-1,0,1,0.1,1,10,4,2,20, 0.0001);

%%%% Sliding
%%%% beam_eq_projection(81,-1,0,1,1,1,10,4,3,20, 0.0001);

%%%% Hinged
%%%% beam_eq_projection(81,-1,0,1,0.1,1,100,4,4,200, 0.00001);

%%%% Stor skillnad i b kräver mindre k. När t.ex. b1=1, b2=100 krävs minst k = 0.00001*h.

function [w,t,x,h,k] = beam_eq_projection_integlobalP(N, x0, xl, xN, T, b1, b2, order, BC, plotspeed, k_ratio, IC, plotornot, timestepperversion)
close all
pause on
%%%Domain%%%
h = (xN-xl)/(N-1);
x1 = x0:h:xl;
x2 = xl:h:xN;
x = [x1 x2];
%%%%%%%%%%%%%%%%%%%%%%%%

%%% Initial condition %%%
if(b1 == b2 && BC == 2 && x0 == 0 && xN == 1) %%%%%%exact solution exists 
   % u0 = @(x) cosh(1.50562*pi.*x) + cos(1.50562*pi.*x) - ((cos(1.50562*pi) -cosh(1.50562*pi))/(sin(1.50562*pi) - sinh(1.50562*pi)))*(sin(1.50562*pi.*x) +sinh(1.50562*pi.*x));
    %u0_t = 0;
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


%%%Boundary conditions Left%%%
if(BC == 1)%clamped
    L_left = [e1; d1_1];
elseif(BC == 2)%free
    L_left = [d1_2; d1_3] ;
elseif(BC == 3)%sliding
    L_left = [d1_1, d1_3]; 
elseif(BC == 4)%hinged
    L_left = [e1; d1_2];
end

%%Boundary Condition Right%%%
if(BC == 1)%clamped
    L_right = [eN; dN_1];
elseif(BC == 2)%free
    L_right = [dN_2; dN_3];
elseif(BC == 3)%sliding
    L_right = [dN_1; dN_3];
elseif(BC == 4)%hinged
    L_right = [eN; dN_2];
end

%%%Inner boundary%%%
L_inner = [eN -e1; dN_1 -d1_1; b1*dN_2 -b2*d1_2; b1*dN_3 -b2*d1_3];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

H_big = [H zeros(N); zeros(N) H];
H_bigI = inv(H_big);
P_left = eye(N) - HI*L_left'*inv(L_left*HI*L_left')*L_left;
P_right = eye(N) - HI*L_right'*inv(L_right*HI*L_right')*L_right;
P_inner = [eye(N) zeros(N); zeros(N) eye(N)] - H_bigI*L_inner'*inv(L_inner*H_bigI*L_inner')*L_inner;
A = (-P_inner*[b1^2*P_left*D4*P_left zeros(N); zeros(N) b2^2*P_right*D4*P_right]*P_inner);


w0 = [u0 v0];
w0_t = u0_t;

if timestepperversion == 1
    [w,k,t] = timestepper(T,h,A,w0,w0_t, k_ratio);
elseif timestepperversion == 2
    [w,k,t] = timestepperv2(T,h,A,w0,w0_t,k_ratio);
end

%%%Exact Solution%%%
if plotornot == 1
    u_exact = @(x,t) real(exp(-1i*(22.3733)*t)*u0(x));
    figure(1);
    for i = 1:plotspeed:length(t)
        if(b1 == b2 && BC == 2 && x0 == 0 && xN == 1)
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