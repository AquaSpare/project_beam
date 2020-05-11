function x = main(BC, a1, a2, ini)

ff = @(x) func(x,a1,a2,BC);
x0 = ones(10,1);
x0(1) = ini;
x0(2) = ini;
x0(3:10) = 10*rand(8,1);
options = optimoptions('fsolve','FunctionTolerance',1e-12);
x = fsolve(ff,x0,options);

ff(x)
