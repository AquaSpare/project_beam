close all

T = 5000;
a1 = 1;
a2 = 4;

N = 41;
order = 4;
k_ratio = 0.19;
            

BC = 1;
x0 = -1;
xl = 0;
xN = 1;

[u_exact_l, u_exact_r] = beam_ana(BC,a1,a2);
error = cell(2,1);
if BC == 1
    [w{1,1},t,x,h,k,error{1,1}] = proj_proj_proj(N, x0, xl, xN, T, a1, a2, order, BC, 100, k_ratio, 1, 0, 3);
    [w{1,2},t,x,h,k,error{2,1}] = proj_projSAT_proj(N, x0, xl, xN, T, a1, a2, order, BC, 100, k_ratio, 1, 0,3);
    
elseif BC == 2
    [w{1,1},t,x,h,k] = SAT_proj_SAT(N, x0, xl, xN, T, a1, a2, order, BC, 100, k_ratio, 1, 0, 1);
    [w{1,2},t,x,h,k] = SAT_projSAT_SAT(N, x0, xl, xN, T, a1, a2, order, BC, 100, k_ratio, 1, 0,1);
    w_exact = zeros(length(t),2*N);
    for i=1:length(t)
    w_exact(i,:) = [u_exact_l(x0:h:xl,(i-1)*k) u_exact_r((xl:h:xN),(i-1)*k)];
    end
end



plot(error{1,1});
hold on
plot(error{2,1});
if BC == 1
    legend('Pure Projection', 'Projection mix');
elseif BC == 2
    legend('SAT-proj','SAT inner mix');
end
hold off



