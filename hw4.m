clear all; close all; clc

% matrix A


m=64; % N value in x and y directions
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

A = spdiags([e1 e1 e5 e2 -4*e1 e3 e4 e1 e1], ...
[-(n-m) -m -m+1 -1 0 1 m-1 m (n-m)],n,n);
% spy(matA) % view the matrix structure

A(1,1) = 2;
A = A/(dn^2);

% matrix B

adds = (1/2/dn)*ones(n, 1);
subs = -1*adds;
B = spdiags([adds adds subs subs], [m -(n-m) (n-m) -m], n, n);


% matrix C
C = spdiags([e5 -e2 e3 -e4], [-(m-1) -1 1 (m-1)], n, n);
C = (1/2/dn).*C;

%%

k = 64;
x = linspace(-10,10,k+1);
y = linspace(-10,10,k+1);

xspan = x(1:k);
yspan = y(1:k);
[X,Y] = meshgrid(xspan,yspan);
tspan = 0:0.5:4;
winit = exp(-X.^2 - (1/20).*(Y.^2));

% surf(winit);
winit = reshape(winit, [k^2 1]);


% [t1,w1] = ode45(@(t1,w1) rhsvs1(t1,w1,A,B,C), tspan, winit);
% 
% [t2,w2] = ode45(@(t2,w2) rhsvs2(t2,w2,A,B,C), tspan, winit);
% 
% [t3,w3] = ode45(@(t3,w3) rhsvs3(t3,w3,A,B,C), tspan, winit);


% temp1 = w1(end,:);
% temp1 = reshape(temp1, 64,64);
% pcolor(temp1);
% 
% temp2 = w2(end,:);
% temp2 = reshape(temp2, 64,64);
% pcolor(temp2);
% 
% temp3 = w3(end,:);
% temp3 = reshape(temp3, 64,64);
% pcolor(temp3);


% save('A1.dat','w1','-ascii');
% save('A2.dat','w2','-ascii');
% save('A3.dat','w3','-ascii');




%%
% Two oppositely “charged” Gaussian vortices next to each other, 
% i.e. one with positive amplitude, the other with negative amplitude.

max_time = 60;
dt = 0.5;
close all;
tspan = 0:dt:max_time;
winit = 0.5*exp(-(X).^2-(Y+1).^2)-0.5*exp(-(X).^2-(Y-1).^2);
winit = reshape(winit, [k^2 1]);
[t,w] = ode45(@(t,w) rhsvs3(t,w,A,B,C), tspan, winit);

loops = size(w,1);
F(loops) = struct('cdata',[],'colormap',[]);
for i = 1:loops
    test = w(i,:);
    test = reshape(test,64,64);
    pcolor(test);
    drawnow;
    F(i) = getframe;
end
% movie(F, 2);
    
v = VideoWriter('m1.avi');
open(v);
writeVideo(v, F);
close(v);

%%
% Two same “charged” Gaussian vortices next to each other.

close all;
tspan = 0:dt:max_time;
winit = 1*exp(-(X).^2-(Y+2).^2)+1*exp(-(X).^2-(Y-2).^2);
winit = reshape(winit, [k^2 1]);
[t,w] = ode45(@(t,w) rhsvs3(t,w,A,B,C), tspan, winit);

loops = size(w,1);
F(loops) = struct('cdata',[],'colormap',[]);
for i = 1:loops
    test = w(i,:);
    test = reshape(test,64,64);
    pcolor(test);
    drawnow;
    F(i) = getframe;
end
% movie(F, 2);
    
v = VideoWriter('m2.avi');
open(v);
writeVideo(v, F);
close(v);

%% 
% Two pairs of oppositely “charged” vortices 
% which can be made to collide with each other.


close all;
tspan = 0:dt:max_time;
winit = 0.8*exp(-(X).^2-(Y+0.5).^2)-1.15*exp(-(X).^2-(Y).^2);
winit = reshape(winit, [k^2 1]);
[t,w] = ode45(@(t,w) rhsvs3(t,w,A,B,C), tspan, winit);

loops = size(w,1);
F(loops) = struct('cdata',[],'colormap',[]);
for i = 1:loops
    test = w(i,:);
    test = reshape(test,64,64);
    pcolor(test);
    drawnow;
    F(i) = getframe;
end
% movie(F, 2);
    
v = VideoWriter('m3.avi');
open(v);
writeVideo(v, F);
close(v);


%%
% A random assortment (in position, strength, charge, ellipticity, etc.) of vortices on the periodic
% domain. Try 10-15 vortices and watch what happens.


close all;
tspan = 0:dt:max_time;
winit = 1*exp(-(X+1).^2-(Y+1).^2)...
    -0.7*exp(-(X).^2-(Y).^2)...
    +1*exp(-(X-1).^2-(Y-1).^2)...
    -1*exp(-(X+1).^2-(Y-1).^2)...
    +0.8*exp(-(X-1).^2-(Y+1).^2)...
    -1*exp(-1/5*(X-4).^2-1/10*(Y+4).^2)...
    +1.1*exp(-1/10*(X-4).^2-1/20*(Y-4).^2)
    -1*exp(-(X+6).^2-(Y).^2)...
    +1*exp(-(X-6).^2-(Y-2).^2)...
    -1.2*exp(-1/2*(X).^2-(Y).^2)...
    +1*exp(-(X-8).^2-(Y-8).^2)...
    +0.5*exp(-1/4*(X+8).^2-1/8*(Y+8).^2)...
    +1*exp(-(X+10).^2-(Y-12).^2);

winit = reshape(winit, [k^2 1]);
[t,w] = ode45(@(t,w) rhsvs3(t,w,A,B,C), tspan, winit);

loops = size(w,1);
F(loops) = struct('cdata',[],'colormap',[]);
for i = 1:loops
    test = w(i,:);
    test = reshape(test,64,64);
    pcolor(test);
    drawnow;
    F(i) = getframe;
end
% movie(F, 2);
    
v = VideoWriter('m_random.avi');
open(v);
writeVideo(v, F);
close(v);

%%

close all;
tspan = 0:dt:max_time;
winit = 1.15*exp(-(X+1).^2-(Y+1).^2)...
    -0.85*exp(-1/5*(X).^2-1/5*(Y).^2)...
    +1*exp(-(X-1).^2-(Y-1).^2);

winit = reshape(winit, [k^2 1]);
[t,w] = ode45(@(t,w) rhsvs3(t,w,A,B,C), tspan, winit);

% test = w(end,:);
% test = reshape(test,64,64);
% figure
% surf(test);
% figure
% pcolor(test);

loops = size(w,1);
F(loops) = struct('cdata',[],'colormap',[]);
for i = 1:loops
    test = w(i,:);
    test = reshape(test,64,64);
    pcolor(test);
    drawnow;
    F(i) = getframe;
end
% movie(F, 2);
    
v = VideoWriter('my_movie.avi');
open(v);
writeVideo(v, F);
close(v);

    
    
    
    
    
    
    
    
    
    
    
    
    
    