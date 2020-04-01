function [D4, H, HI, M4, e1, eN, d1_1, dN_1, d1_2, dN_2, d1_3, dN_3] = SBP2_D4(m, h)

H=diag(ones(m,1),0);H(1,1)=1/2;H(m,m)=1/2;


H=H*h;
HI=inv(H);
e1=zeros(1,m);e1(1)=1;
eN=zeros(1,m);eN(m)=1;
m2=1;m1=-4;m0=6;
M4=m2*(diag(ones(m-2,1),2)+diag(ones(m-2,1),-2))+m1*(diag(ones(m-1,1),1)+diag(ones(m-1,1),-1))+m0*diag(ones(m,1),0);

%M4=(-1/6*(diag(ones(m-3,1),3)+diag(ones(m-3,1),-3) ) + 2*(diag(ones(m-2,1),2)+diag(ones(m-2,1),-2)) -13/2*(diag(ones(m-1,1),1)+diag(ones(m-1,1),-1)) + 28/3*diag(ones(m,1),0));

M4_U=[0.13e2 / 0.10e2 -0.12e2 / 0.5e1 0.9e1 / 0.10e2 0.1e1 / 0.5e1; -0.12e2 / 0.5e1 0.26e2 / 0.5e1 -0.16e2 / 0.5e1 0.2e1 / 0.5e1; 0.9e1 / 0.10e2 -0.16e2 / 0.5e1 0.47e2 / 0.10e2 -0.17e2 / 0.5e1; 0.1e1 / 0.5e1 0.2e1 / 0.5e1 -0.17e2 / 0.5e1 0.29e2 / 0.5e1;];


M4(1:4,1:4)=M4_U;

M4(m-3:m,m-3:m)=flipud( fliplr( M4_U ) );
M4=M4/h^3;

d1_U=[-3/2 2 -1/2]/h;
d1_1=zeros(1,m);
d1_1(1:3)=d1_U;
dN_1=zeros(1,m);
dN_1(m-2:m)=fliplr(-d1_U);

d2_U=[1 -2 1;]/h^2;
d1_2=zeros(1,m);
d1_2(1:3)=d2_U;
dN_2=zeros(1,m);
dN_2(m-2:m)=fliplr(d2_U);

d3_U=[-1 3 -3 1;]/h^3;
d1_3=zeros(1,m);
d1_3(1:4)=d3_U;
dN_3=zeros(1,m);
dN_3(m-3:m)=fliplr(-d3_U);

D4=HI*(M4-e1'*d1_3+eN'*dN_3  + d1_1'*d1_2-dN_1'*dN_2);