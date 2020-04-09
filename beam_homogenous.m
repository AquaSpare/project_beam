
pause on
clear all
close all

%mesh and settings
N = 201;
x1 = 0;
x2 = 1;
h = (x2-x1)/(N-1);
x = x1:h:x2;
t0 = 0;
T = 0.01;

%constant
b = 1;
xp = 0.5;
r0 = 0.05;

%initial condition
u0 = @(x) cosh(1.50562*pi.*x) + cos(1.50562*pi.*x) - ((cos(1.50562*pi) -cosh(1.50562*pi))/(sin(1.50562*pi) - sinh(1.50562*pi)))*(sin(1.50562*pi.*x) +sinh(1.50562*pi.*x));
w0_t = 0;
w0 = u0(x);


[D4, H, HI, M, e1, eN, d1_1, dN_1, d1_2, dN_2, d1_3, dN_3] = SBP4_D4(N, h);
HI = inv(H);

%------Clamped-------
%L = [e1;eN;d1_1;dN_1];
%------Hinged--------
%L = [e1;eN;d1_2;dN_2];
%------Sliding-------
%L = [d1_1;dN_1;d1_3;dN_3];
%-------Free---------
L = [d1_2;dN_2;d1_3;dN_3];

P = eye(N)-HI*L'*inv(L*HI*L')*L;
A = -b.*P*D4*P;

[w, k, t] = timestepper(t0, T, h, A, w0, w0_t, 0.0001);

for i = 1:length(t)
    plot(x,w(:,i))
    axis([0 1 -1 5])
    pause(0.00001)
end

