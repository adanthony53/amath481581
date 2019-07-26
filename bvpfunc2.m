function dydt = bvpfunc2(x,y,k,gamma,eps)
f1 = y(2);
f2 = (gamma*(abs(y(1))^2) + k*x^2-eps)*y(1);
dydt = [f1;f2];
end