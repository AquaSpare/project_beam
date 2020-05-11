function [e1,e2,steps] = homo_konvergens(iter,BC,order,a)

x0 = 0;
xN = 1;

N = 14;
T = 0.05;
kratio = 0.0001;

e1 = zeros(1,iter);
e2 = zeros(1,iter);

steps = zeros(1,iter);

w_ana = homo_beam_ana_2(a, BC);

for i = 1:iter
    [wsat, tsat, ~, ~, ~, ~] = SAT_enBalk((N*i)+1,x0,xN,T, a, order, BC, 1000, kratio, 2, 0, 2);
    [wproj,tproj,x,h,~,~] = beam_homogenous((N*i)+1,x0,xN,T,a,order,BC,1000,kratio,2,0, 2);
    
    steps(i) = h;
    
    e1(i) = sqrt(1/(N*i))*norm(w_ana(x,tsat(end))-wsat(:,end)',2);
    e2(i) = sqrt(1/(N*i))*norm(w_ana(x,tproj(end))-wproj(:,end)',2);
end

plot(steps,e1,'r*')
hold on
plot(steps,e2,'b*')
legend('SAT','Projection')
