function [D2, H, HI, e1, eN, d1, dN] = SBP4(N, h)

e1 = zeros(1,N);
eN = zeros(1,N);
d1 = zeros(1,N);
dN = zeros(1,N);
e1(1,1) = 1;   
eN(1,end) = 1;
d1(1,1:4) = 1/h*[11/6 -3 3/2 -1/3];
dN(1,(end-3):end) = 1/h*[1/3 -3/2 3 -11/6];

H = diag(ones(1,N));
H(1,1) = 17/48; H(end,end) = 17/48;
H(2,2) = 59/48; H(end-1,end-1) = 59/48;
H(3,3) = 43/48; H(end-2,end-2) = 43/48;
H(4,4) = 49/48; H(end-3,end-3) = 49/48;
H = h * H;
HI = inv(H);

M = (zeros(N)  + diag(5/2*ones(N ,1)) + diag(-4/3*ones(N-1,1), 1) + diag(-4/3*ones(N-1,1), -1) + diag(1/12*ones(N-2,1), 2) + diag(1/12*ones(N-2,1), -2));  
M_U1 = [9/8 -59/48 1/12 1/48]; M_U2 = [-59/48 59/24 -59/48 0]; M_U3 = [1/12 -59/48 55/24 -59/48]; M_U4 = [1/48 0 -59/48 59/24];
M(1,1:4) = M_U1; M(2,1:4) = M_U2; M(3,1:4) = M_U3; M(4,1:4) = M_U4;
M(end,(end-3):end) = fliplr(M_U1); M(end-1,(end-3):end) = fliplr(M_U2); M(end-2,(end-3):end) = fliplr(M_U3); M(end-3,(end-3):end) = fliplr(M_U4);
M = 1/h * M;

D2 = H\(-M-e1'*d1+eN'*dN);