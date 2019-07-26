clear all; close all; clc

% y = @(t) pi*exp(3*(cos(t)-1))/sqrt(2);
% f = @(t, y) -3 * y * sin(t);
%
% % dt = [2^-10 2^-8];
%
% dt = [2^-2 2^-3 2^-4 2^-5 2^-6 2^-7 2^-8];
% E1 = [];
% E2 = [];
%
% for delta_t = dt
%     t = [0:delta_t:5];
%     y1 = [pi/sqrt(2)];
%     y2 = [pi/sqrt(2)];
%     y_true = y(t);
%
%     iter = 2;
%     for i = t(2:end)
%
%         y1(iter) = y1(iter-1) + delta_t * f(i, y1(iter-1));
%
%         y2(iter) = y2(iter-1) + delta_t/2 * (...
%         f(i, y2(iter-1)) + ...
%         f(i+delta_t, y2(iter-1)+delta_t*f(i, y2(iter-1)))...
%         );
%
%         iter = iter + 1;
%     end
% %     hold on
% %     plot(t, y_true, 'r');
% %     plot(t, y1, 'b');
% %     plot(t, y2, 'g');
%
%     E1 = [E1 mean(abs(y_true - y1))];
%     E2 = [E2 mean(abs(y_true - y2))];
% end
% % plot(log(dt), log(E1));
% % p1 = polyfit(log(dt), log(E1), 1);
% % p2 = polyfit(log(dt), log(E2), 1);
% E1
%
% a1 = transpose(y1);
% a2 = E1;
% a3 = p1(1);
% % save('A1.dat','a1','-ascii');
% % save('A2.dat','a2','-ascii');
% % save('A3.dat','a3','-ascii');
% a4 = transpose(y2);
% a5 = E2;
% a6 = p2(1);
% % save('A4.dat','a4','-ascii');
% % save('A5.dat','a5','-ascii');
% % save('A6.dat','a6','-ascii');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


dt = [2^-2 2^-3 2^-4 2^-5 2^-6 2^-7 2^-8];

f = @(t, y) -3 * y * sin(t);
y = @(t) pi*exp(3*(cos(t)-1))/sqrt(2);
E1 = [];
E2 = [];

for k = [1:length(dt)]
    h = dt(k);
    t = [0:h:5];
    y_true = y(t);
    y1 = zeros(size(y_true)); y1(1) = pi/sqrt(2);
    y2 = zeros(size(y_true)); y2(1) = pi/sqrt(2);
    
    for i = 2:length(t)
        %ti = t(i-1);
        
        y1(i) = y1(i-1) + h*f(t(i-1), y1(i-1));
        
        y2(i) = y2(i-1) + h/2*(f(t(i-1), y2(i-1)) + ...
            f(t(i-1)+h, y2(i-1)+h*f(t(i-1), y2(i-1))));
    end
    E1 = [E1 mean(abs(y_true - y1))];
    E2 = [E2 mean(abs(y_true - y2))];
    
    % close all;
    % hold on
    % plot(t, y_true, 'r');
    % plot(t, y1, 'b');
    % plot(t, y2, 'g');
    
end

p1 = polyfit(log(dt), log(E1), 1);
p2 = polyfit(log(dt), log(E2), 1);

a1 = transpose(y1);
a2 = E1;
a3 = p1(1);
save('A1.dat','a1','-ascii');
save('A2.dat','a2','-ascii');
save('A3.dat','a3','-ascii');
a4 = transpose(y2);
a5 = E2;
a6 = p2(1);
save('A4.dat','a4','-ascii');
save('A5.dat','a5','-ascii');
save('A6.dat','a6','-ascii');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc

tspan = [0:0.5:32];
yinit = [sqrt(3); 1];
eps = 0.1;
[T1, Y1] = ode45(@(t,y) ode1(t,y,eps), tspan, yinit);
eps = 1;
[T2, Y2] = ode45(@(t,y) ode1(t,y,eps), tspan, yinit);
eps = 20;
[T3, Y3] = ode45(@(t,y) ode1(t,y,eps), tspan, yinit);

a7 = [Y1(:,1) Y2(:,1) Y3(:,1)];
save('A7.dat','a7','-ascii');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear all; close all; clc


tspan = [0 32];
yinit = [2; pi*pi];
eps = 1;
tol_val = [1e-4 1e-5 1e-6 1e-7 1e-8 1e-9 1e-10];

step_size_45 = [];
step_size_23 = [];
step_size_113 = [];

for tol = tol_val
    options = odeset('AbsTol',tol,'RelTol',tol);
    
    [T1, Y1] = ode45(@(t,y) ode1(t,y,eps), tspan, yinit, options);
    [T2, Y2] = ode23(@(t,y) ode1(t,y,eps), tspan, yinit, options);
    [T3, Y3] = ode113(@(t,y) ode1(t,y,eps), tspan, yinit, options);
    
    step_size_45 = [step_size_45 mean(diff(T1))];
    step_size_23 = [step_size_23 mean(diff(T2))];
    step_size_113 = [step_size_113 mean(diff(T3))];
end

% figure
% hold on
% plot(log(step_size_45), log(tol_val), '-ro');
% plot(log(step_size_23), log(tol_val), '-bo');
% plot(log(step_size_113), log(tol_val), '-go');

pode45 = polyfit(log(step_size_45), log(tol_val), 1);
pode23 = polyfit(log(step_size_23), log(tol_val), 1);
pode113 = polyfit(log(step_size_113), log(tol_val), 1);

a8 = pode45(1);
a9 = pode23(1);
a10 = pode113(1);

save('A8.dat','a8','-ascii');
save('A9.dat','a9','-ascii');
save('A10.dat','a10','-ascii');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc

a1 = 0.05;
a2 = 0.25;
b = 0.01;
c = 0.01;
I = 0.1;

tspan = [0:0.5:100];
init = [0.1; 0.1; 0; 0];

%a11 (0,0)
[t11,y11] = ode15s(@(t,y) fitzhugh_coupling(t,y,a1,a2,b,c,I,0,0), ...
    tspan, init);

%a12 (0,0.2)
[t12,y12] = ode15s(@(t,y) fitzhugh_coupling(t,y,a1,a2,b,c,I,0,0.2), ...
    tspan, init);

%a13 (-0.1,0.2)
[t13,y13] = ode15s(@(t,y) fitzhugh_coupling(t,y,a1,a2,b,c,I,-0.1,0.2), ...
    tspan, init);

%a14 (-0.3,0.2)
[t14,y14] = ode15s(@(t,y) fitzhugh_coupling(t,y,a1,a2,b,c,I,-0.3,0.2), ...
    tspan, init);

%a15 (-0.5,0.2)
[t15,y15] = ode15s(@(t,y) fitzhugh_coupling(t,y,a1,a2,b,c,I,-0.5,0.2), ...
    tspan, init);


save('A11.dat','y11','-ascii');
save('A12.dat','y12','-ascii');
save('A13.dat','y13','-ascii');
save('A14.dat','y14','-ascii');
save('A15.dat','y15','-ascii');









