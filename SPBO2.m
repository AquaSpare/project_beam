
function [H,M,e1,em,d1,dm,D2] = SPBO2(N,h)

H = diag(ones(N,1));
H(1,1) = 1/2;
H(end,end) = 1/2;
H = h.*H;

M = diag(2.*ones(N,1)) + diag(-1.*ones(N-1,1),1) + diag(-1.*ones(N-1,1),-1);
M(1,1) = 1;
M(end,end) = 1;
M = h.*M;

e1 = zeros(1,N);
e1(1,1) = 1;
em = zeros(1,N);
em(1,end) = 1;

d1 = zeros(1,N);
d1(1,1:3) = [-3/2 2 1/2];
d1 = d1.*(1/h);
dm = zeros(1,N);
dm(1,(end-2):end) = [1/2 -2 3/2];
dm = dm.*(1/h);

D2 = H\(-M-e1*d1+em*dm);
end
