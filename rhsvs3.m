function dwdt = rhsvs3(t3,w3,A,B,C)
k = 64;
wmatrix = reshape(w3, k, k);
kx = (2*pi/20)*[0:(k-1)/2 -k/2:-1];
ky = (2*pi/20)*[0:(k-1)/2 -k/2:-1];
kx(1) = 10e-6;
ky(1) = 10e-6;
[KX,KY] = meshgrid(kx,ky);

phimatrix = real(ifft2(-fft2(wmatrix)./(KX.^2+KY.^2)));
phi = reshape(phimatrix, k^2, 1);

dwdt = -(B*phi).*(C*w3) + (C*phi).*(B*w3) + 0.001*A*w3;

end
