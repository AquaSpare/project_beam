close all;
clear all;
N = 41; x0 = -1; xl =0; xN = 1; a1 = 1; a2 = 100; order = 4; BC = 2; k_ratio = 0.00095;

h = (xN-xl)/(N-1);
x1 = x0:h:xl;
x2 = xl:h:xN;
x = [x1 x2];
k = h*k_ratio;
%%%%%%%%%%%%%%%%%%%%%%%%


%%%SBP operators%%%
if(order == 2)
    [D4, H, HI, M, e1, eN, d1_1, dN_1, d1_2, dN_2, d1_3, dN_3] = SBP2_D4(N, h);
elseif(order == 4)
    [D4, H, HI, M, e1, eN, d1_1, dN_1, d1_2, dN_2, d1_3, dN_3] = SBP4_D4(N, h);
elseif(order == 6)
    [D4, H, HI, M, e1, eN, d1_1, dN_1, d1_2, dN_2, d1_3, dN_3] = SBP6_D4(N, h);
end


tau5 = linspace(0.99,100,100);
tau6 = 1;
tau3 = 1 + (tau5-1)*a2/a1;
tau4 = 1 + (tau6-1)*a2/a1;
maxeig = zeros(length(tau5),1);

L = [eN -e1; dN_1 -d1_1];
H = [H zeros(N); zeros(N) H];
HI2 = inv(H);
P = [eye(N) zeros(N); zeros(N) eye(N)] - HI2*L'*inv(L*HI2*L')*L;

for i=1:length(tau5)
   %for j=1:length(tau6)
        l_u = -a1*D4 - a1*HI*(-d1_1'*d1_2 + e1'*d1_3 - tau3(i)*eN'*dN_3 + tau4*dN_1'*dN_2);
        r_l = -a2*D4 - a2*HI*(-tau6*d1_1'*d1_2 + tau5(i)*e1'*d1_3 - eN'*dN_3 + dN_1'*dN_2);
        A = P*[l_u zeros(N); zeros(N) r_l]*P;
        maxeig(i) = max(abs(eig(k^2*A+2*eye(size(A)))));
   %end
end
plot(tau5,maxeig);
xlabel('Tau5');
ylabel('Tau6');
zlabel('Maxeig');
