function dwdt = rhsvs2(t,w2,A,B,C)
[L,U,P] = lu(A);
y = L\(P*w2);
phi = U\y;
dwdt = -(B*phi).*(C*w2) + (C*phi).*(B*w2) + 0.001*A*w2;
end
