%function [w,t,x,h,k] = beam_eq_projection_dav(N,x0,xl,xN,t0,T,b1,b2,order,BC)

close all

[u,v] = beam_ana();

[w,t,x,h,k] =  beam_eq_projection(81,-1,0,1,0.1,1,4,4,2,20, 0.0001);
N = 81
x1 = 0:h:N;
x2 = N:h:2*N;
w_ana = [u(x1,t) v(x2,t)];

norm(w-w_ana,2)