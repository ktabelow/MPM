function [t,iter] = Armijo(fctn,dx,xk,f,df)

c1 = 1e-3;
maxIter = 100;
iter = 1;

t = 1;
beta = 0.7;
while iter <= maxIter
    xt = xk + t*dx;
    ft = fctn(xt);
    
    if (ft <= f + c1*df'*dx)
        return
    end
    t = t*beta;
    iter = iter + 1;
end
    
    