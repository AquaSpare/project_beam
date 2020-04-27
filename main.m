function x = main(BC, a1, a2, ini)

ff = @(x) func(x,a1,a2,BC);
x0 = ones(10,1);
x0(1) = ini;
x0(2) = ini;
options = optimoptions('fsolve','FunctionTolerance',1e-15);
x = fsolve(ff,x0,options);

ff(x)
