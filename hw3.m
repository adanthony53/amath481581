clear all; close all; clc
clear all; close all; % clear all variables and figures

m=4; % N value in x and y directions
n=m*m; % total size of matrix
dn = 20/m;

e0=zeros(n,1); % vector of zeros
e1=ones(n,1); % vector of ones

e2=e1; % copy the one vector
e4=e0; % copy the zero vector

for j=1:m
    e2(m*j)=0; % overwrite every m^th value with zero
    e4(m*j)=1; % overwirte every m^th value with one
end

e3(2:n,1)=e2(1:n-1,1); 
e3(1,1)=e2(n,1); % shift to correct
e5(2:n,1)=e4(1:n-1,1); 
e5(1,1)=e4(n,1); % positions
% place diagonal elements

matA = spdiags([e1 e1 e5 e2 -4*e1 e3 e4 e1 e1], ...
[-(n-m) -m -m+1 -1 0 1 m-1 m (n-m)],n,n);
% spy(matA) % view the matrix structure
A(1,1) = 2;
A = A./(dn^2);

%% matrix B


adds = (1/2/dn)*ones(n, 1);
subs = -1*adds;
B = full(spdiags([adds adds subs subs], [m -(n-m) (n-m) -m], n, n));



%% matrix C
e0 = zeros(n, 1);
e1 = ones(n, 1);
e2 = e0;
e3 = e1;

for i=1:n
    if mod(i, m) == 1
        e0(i) = 1;
        e1(i) = 0;
    end
    if mod(i, m) == 0
        e2(i) = 1;
        e3(i) = 0;
    end
end

e2 = -1*e0;
e3 = -1*e3;
C = full(spdiags([e1 e3 e0],[1 -1 -7],n,n));

for i=1:m
    C(1+(i-1)*8, 8*i) = -1;
end
A3 = C*(1/2/dn);


%%

% clear all; close all; clc
fmat = load('Fmat.mat');
permvec = load('permvec.mat');

n = fieldnames(fmat);
m = getfield(fmat,n{1});

n = fieldnames(permvec);
v = getfield(permvec,n{1});

center = m(161:240,161:240);
% uint8(abs(center))


C = mat2cell(center, [20 20 20 20], [20 20 20 20]);
C = reshape(C, 1, 16); %4*4 to 1*16

%%%%%%%%%%%%%%%%%%%%%%
% 7     4     2    13%
%11    16     6     9%
% 3     5     1     8%
%14     15   10    12%
%%%%%%%%%%%%%%%%%%%%%%


D = C;
for i = 1:length(v)
    D(i) = C(v(i));
end
D = reshape(D, 4, 4);
D = cell2mat(D);


% D = cell2mat([C(3,2), C(4,1), C(2,1), C(1,4);...
%               C(3,3), C(4,4), C(2,1), C(1,3);...
%               C(3,1), C(1,2), C(1,1), C(4,2);...
%               C(2,4), C(3,4), C(2,3), C(4,3)]);
          
m(161:240,161:240) = D;

A4 = abs(m);

A5 = abs(ifft2(ifftshift(m)));



save('A1.dat','A1','-ascii');
save('A2.dat','A2','-ascii');
save('A3.dat','A3','-ascii');
save('A4.dat','A4','-ascii');
save('A5.dat','A5','-ascii');




















