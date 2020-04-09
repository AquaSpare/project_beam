pause on
clear all
close all

%mesh and settings
N = 101;
x1 = 0;
x2 = 1;
h = (x2-x1)/(N-1);
x = x1:h:x2;
t0 = 0;
T = 1;

%constant
c = 4;
xp = 0.5;
r0 = 0.05;

%initial data
u0 = @(x) exp(-(xp-x).^2/(r0^2));
w0_t = 0;
w0 = u0(x);


[D2, H, HI, M, e1, eN, d1, dN] = SBP4(N, h);
HI = inv(H);

L = [d1;eN];
P = eye(N)-HI*L'*inv(L*HI*L')*L;
A = (c^2)*P*D2*P;

[w, k, t] = timestepper(t0, T, h, A, w0, w0_t);

for i = 1:length(t)
    plot(x,w(:,i))
    axis([0 1 -1 1])
    pause(0.0000000001);
end

