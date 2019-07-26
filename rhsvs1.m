function dwdt = rhsvs1(t1,w1,A,B,C)
phi = A\w1;
% dwdt = -1*dot(B*phi,C*w1) + dot(C*phi,B*w1) + 0.001.*(A*w1);
dwdt = -(B*phi).*(C*w1) + (C*phi).*(B*w1) + 0.001*A*w1;
end
