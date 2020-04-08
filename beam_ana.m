% Provides analytical solution to the problem
% utt = -a1*uxxxx, -1 < x < 0,
% vtt = -a2*vxxxx, 0 < x < 1.
% where a1 = 1 and a2 = 4.
%
% Boundary conditions: 
% uxx = uxxx = 0, at x = -1,
% vxx = vxxx = 0, at x = 1. 
%
% Coupling conditions:
% u = v, ux = vx, a1*uxx = a2*vxx, a1*uxxx = a2*vxxx, at x = 0.
%
% Initial conditions:
% u = X1, v = X2, ut = 0, vt = 0, at t = 0
%
function [u,v] = beam_ana()
a1 = 1;
a2 = 4;

params = [2.636119548106415;
   1.864018008484463;
   1.069297874524326;
  -0.286268581042760;
  -0.256470338226990;
  -0.144340958704486;
   0.737855821336497;
  -0.354665917695434;
   0.074971714960839;
  -0.254307933502358];
   
w1 = params(1);
w2 = params(2);
b3 = params(3);
b4 = params(4);
b5 = params(5);
b6 = params(6);
c3 = params(7);
c4 = params(8);
c5 = params(9);
c6 = params(10);

T1 = @(t) cos(sqrt(a1)*w1^2*t);
X1 = @(x) b3*cos(w1*x) + b4*sin(w1*x) + b5*cosh(w1*x) + b6*sinh(w1*x);
T2 = @(t) cos(sqrt(a2)*w2^2*t);
X2 = @(x) c3*cos(w2*x) + c4*sin(w2*x) + c5*cosh(w2*x) + c6*sinh(w2*x);

u = @(x,t) T1(t)*X1(x);
v = @(x,t) T2(t)*X2(x);