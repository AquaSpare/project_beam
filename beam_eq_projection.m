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

function [w,t,x,h,k] = beam_eq_projection(N, x0, xl, xN, T, b1, b2, order, BC, plotspeed, k_ratio)
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
    u0 = @(x) cosh(1.50562*pi.*x) + cos(1.50562*pi.*x) - ((cos(1.50562*pi) -cosh(1.50562*pi))/(sin(1.50562*pi) - sinh(1.50562*pi)))*(sin(1.50562*pi.*x) +sinh(1.50562*pi.*x));
    u0_t = 0;
else 
    xr = -1/4;
    r0 = 1/30;
    u0 = @(x) exp(-(xr-x).^2/r0^2);
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
A = (-P*[b1*D4 zeros(N); zeros(N) b2*D4]*P);


w0 = [u0(x1) u0(x2)];
w0_t = u0_t;

[w,k,t] = timestepper(T,h,A,w0,w0_t, k_ratio);

%%%Exact Solution%%%
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