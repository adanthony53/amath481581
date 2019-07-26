% v1    v2       w1      w2 
% y(1)  y(2)     y(3)    y(4)


function dydt = fitzhugh_coupling(t,y,a1,a2,b,c,I,d1,d2)

v1 = y(1);
v2 = y(2);

w1 = y(3);
w2 = y(4);

dv1 = -v1^3 + (1+a1).*v1.^2 - a1*v1 - w1 + I + d1*v2;
dw1 = b*v1 - c*w1;

dv2 = -v2^3 + (1+a2).*v2.^2 - a2*v2 - w2 + I + d2*v1;
dw2 = b*v2 - c*w2;

dydt = [dv1;dv2;dw1;dw2];
end

