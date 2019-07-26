function dydt = bvpfunc(x,y,k,eps)
f1 = y(2);
f2 = (k*x^2-eps)*y(1);
dydt = [f1;f2];
end