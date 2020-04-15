function [w,t,x,h,k] = beam_test_proj(N,x0,xl,xN,T,b1,b2,order,BC,plotspeed,kratio,plotornot)
close all
pause on

h = (xl-x0)/(N-1);
x1 = x0:h:xl;
x2 = xl:h:xN;
x = [x1 x2];

%initial conditions
if b1 == b2 && BC == 4 && x0 == 0 && xN == 1
    u0 = @(x) cosh(1.50562*pi.*x) + cos(1.50562*pi.*x) - ((cos(1.50562*pi) -cosh(1.50562*pi))/(sin(1.50562*pi) - sinh(1.50562*pi)))*(sin(1.50562*pi.*x) +sinh(1.50562*pi.*x));
    w0_t = 0;
    w0 = u0(x);
    special_case = 1;
else
    xp = -1/4;
    r0 = 1/20;
    u0 = @(x) exp(-(xp-x).^2/r0^2);
    w0_t = 0;
    special_case = 0;
end

if(order == 2)
    [D4, H, HI, M, e1, eN, d1_1, dN_1, d1_2, dN_2, d1_3, dN_3] = SBP2_D4(N, h);
elseif(order == 4)
    [D4, H, HI, M, e1, eN, d1_1, dN_1, d1_2, dN_2, d1_3, dN_3] = SBP4_D4(N, h);
elseif(order == 6)
    [D4, H, HI, M, e1, eN, d1_1, dN_1, d1_2, dN_2, d1_3, dN_3] = SBP6_D4(N, h);
end


if BC == 1
    %------Clamped-------
    L = [eN -e1; dN_1 -d1_1; b1*dN_2 -b2*d1_2; b1*dN_3 -b2*d1_3; e1 zeros(1,N); d1_1 zeros(1,N); zeros(1,N) eN; zeros(1,N) dN_1];
elseif BC == 2
    %------Hinged--------
    L = [eN -e1; dN_1 -d1_1; b1*dN_2 -b2*d1_2; b1*dN_3 -b2*d1_3;e1 zeros(1,N); d1_2 zeros(1,N); zeros(1,N) eN; zeros(1,N) dN_2];
elseif BC == 3
    %------Sliding-------
    L = [eN -e1; dN_1 -d1_1; b1*dN_2 -b2*d1_2; b1*dN_3 -b2*d1_3;d1_1 zeros(1,N); d1_3 zeros(1,N); zeros(1,N) dN_1; zeros(1,N) dN_3];
elseif BC == 4
    %-------Free---------
    L = [eN -e1; dN_1 -d1_1; b1*dN_2 -b2*d1_2; b1*dN_3 -b2*d1_3;d1_2 zeros(1,N); d1_3 zeros(1,N); zeros(1,N) dN_2; zeros(1,N) dN_3];
end

H = [H zeros(N); zeros(N) H];
HI = inv(H);

P = [eye(N) zeros(N); zeros(N) eye(N)] - HI*L'*inv(L*HI*L')*L;
A = P*[-b1*D4 zeros(N); zeros(N) -b2*D4]*P;

w0 = [u0(x1) u0(x2)];
w0_t = 0;

[w,k,t] = timestepper(T,h,A,w0,w0_t, kratio);

if plotornot == 1
    for i = 1:plotspeed:length(t)
        plot(x1, w(1:N,i), '*b', x2, w(N+1:end,i), '*r');
        axis([x0 xN -1 1]);
        hold off
        pause(0.00000001);
    end
end

end