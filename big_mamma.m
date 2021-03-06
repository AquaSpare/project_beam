close all
homoornot = 0;

b1 = 1;
b2 = 4;
order = 4;
BC = 1;
N = 41;
resolution = 1000;
eigs = zeros(5,resolution);

if homoornot == 0
    ratios = linspace(0.194,0.195,resolution);
    for i = 1:3
        for j = 1:length(ratios)  
            eigs(i,j) = eigmax(i,N,ratios(j),BC,order,b1,b2); 
        end
    end
    plot(ratios,eigs(1,:))
    hold on
    plot(ratios,eigs(2,:))
    hold on
    plot(ratios,eigs(3,:))
    axis([0.194 0.195 1.8 2.2])
    legend('Pure projection','proj\_projSAT\_proj','SAT\_projSAT\_SAT')
    xlabel('ratio')
    ylabel('spectral radius')
    title('Spectral radius of diefferent discretization ratios')
    hold off
elseif homoornot == 1
    ratios = linspace(0.3,0.45,1000);
    for i = 4:5
        for j = 1:length(ratios)  
            eigs(i,j) = eigmax(i,N,ratios(j),BC,order,b1,b2); 
        end
    end
    plot(ratios,eigs(4,:))
    hold on
    plot(ratios,eigs(5,:))
    xlabel('ratio')
    ylabel('spectral radius')
    title('Spectral radius for different discretization ratios')
    legend('SAT','projection')
    hold off
end



