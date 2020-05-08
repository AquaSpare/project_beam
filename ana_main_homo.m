function x = ana_main_homo(BC, ini)

ff = @(x) ana_func_homo(x, BC);
x0 = ones(5,1);
x0(1) = ini;
for i=2:5
    x0(i) = 10*rand(1);
end
options = optimoptions('fsolve','FunctionTolerance',1e-15);
x = fsolve(ff,x0,options);

ff(x)
