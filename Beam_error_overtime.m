close all

T = 5;
a1 = 1;
a2 = 4;

N = 161;
order = 4;
k_ratio = 0.0006;
            

BC = 1;
x0 = -1;
xl = 0;
xN = 1;

[u_exact_l, u_exact_r] = beam_ana(BC,a1,a2);

if BC == 1
    [w{1,1},t,x,h,k] = proj_proj_proj(N, x0, xl, xN, T, a1, a2, order, BC, 100, k_ratio, 1, 0, 1);
    [w{1,2},t,x,h,k] = proj_projSAT_proj(N, x0, xl, xN, T, a1, a2, order, BC, 100, k_ratio, 1, 0,1);
    w_exact = zeros(length(t),2*N);
    for i=1:length(t)
    w_exact(i,:) = [u_exact_l(x0:h:xl,(i-1)*k) u_exact_r((xl:h:xN),(i-1)*k)];
    end
    
elseif BC == 2
    [w{1,1},t,x,h,k] = SAT_proj_SAT(N, x0, xl, xN, T, a1, a2, order, BC, 100, k_ratio, 1, 0, 1);
    [w{1,2},t,x,h,k] = SAT_projSAT_SAT(N, x0, xl, xN, T, a1, a2, order, BC, 100, k_ratio, 1, 0,1);
    w_exact = zeros(length(t),2*N);
    for i=1:length(t)
    w_exact(i,:) = [u_exact_l(x0:h:xl,(i-1)*k) u_exact_r((xl:h:xN),(i-1)*k)];
    end
end

error = zeros(2,length(t));
for i = 1:length(t)
error(1,i) = sqrt(h)*norm(w_exact(i,:)' - w{1,1}(:,i),2);
error(2,i) = sqrt(h)*norm(w_exact(i,:)' - w{1,2}(:,i),2);
end

plot(t,error(1,:));
hold on
plot(t,error(2,:));
if BC == 1
    legend('Pure Projection', 'Projection mix');
elseif BC == 2
    legend('SAT-proj','SAT inner mix');
end



