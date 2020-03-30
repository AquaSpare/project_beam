
h = 0.1;
k = 1;
m = 1/h + 1;

H = diag(ones(m));
H(1,1) = 1/2;
H(end,end) = 1/2;
H = h.*H;

