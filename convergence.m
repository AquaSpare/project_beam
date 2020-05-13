%-------------FUNKAR INTE MED SUBPLOTTARNA ----------------------


%convergence
%program för att generera homogen konvergens plot
close all

iter = 2;
order = 2;
a = 1;

e_sat = zeros(4,iter);
e_proj = zeros(4,iter);

figure(1)
for i = 1:4
    [e_sat(i,:),e_proj(i,:),steps] = homo_konvergens(iter,i,order,a,0);
    
    subplot(2,2,i)
    plot(steps,e_sat(i,:),'rx')
    hold on
    plot(steps,e_proj(i,:),'b*')
    legend('SAT','Projection')
    hold off
    
    if i == 1
        title('Homogenous beam, clamped boundary')
    elseif i == 2
        title('Homogenous beam, free boundary')
    elseif i == 3
        title('Homogenous beam, hinged boundary')
    elseif i == 4
        title('Homogenous beam, sliding boundary')
    else
        return
    end
end