clear all; close all; clc

format long 
%%

eps_start = 0;
eps = eps_start;
% deps = 0.001;

xspan = [-4:0.1:4];
tol = 1e-4;
k = 1;

eig_val1 = [];
eig_vec1 = [];

for mode = 1:5
    eps = eps_start;
    deps = 0.1;
    
    for i = 1:1000
        [x,y] = ode45(@(x,y) bvpfunc(x,y,k,eps), xspan, [1;sqrt(16-eps)]);
        temp = y(end,2)+sqrt(16-eps)*y(end,1);
        
        if abs(temp) < tol 
            eig_val1 = [eig_val1; eps];
            y_norm = y(:,1)./sqrt(trapz(xspan, y(:,1).^2));
            eig_vec1 = [eig_vec1 abs(y_norm)];
            break;
        end
        
        % y = y(:,1) / norm(y(:,1));
        % c = trapz(xspan, y.^2);
        % y_norm = y/sqrt(c);
        
        if (-1)^(mode+1)*temp > 0
            eps = eps + deps;
        else
            eps = eps - deps/2;
            deps = deps/2;
        end
    end
    eps_start = eps + 0.1;
end

A1 = eig_vec1(:,1);
A2 = eig_vec1(:,2);
A3 = eig_vec1(:,3);
A4 = eig_vec1(:,4);
A5 = eig_vec1(:,5);
A6 = eig_val1;

save('A1.dat','A1','-ascii');
save('A2.dat','A2','-ascii');
save('A3.dat','A3','-ascii');
save('A4.dat','A4','-ascii');
save('A5.dat','A5','-ascii');
save('A6.dat','A6','-ascii');

%%
clear all; close all; clc

dx = 0.1;
xspan = [-4:dx:4];
len = length(xspan);

A = zeros(79,79); % 79*79
A(1,1) = 2/3+dx*dx*xspan(2)^2; 
A(1,2) = -2/3;
for i = 2:78
    A(i, i-1) = -1;
    A(i, i) = 2+dx*dx*xspan(i+1)^2;
    A(i, i+1) = -1;
end
A(79,78) = -2/3;
A(79,79) = 2/3 + dx*dx*xspan(len-1)^2;

[V, D] = eig(A);
[eig_val2, index] = sort(diag(D)/dx/dx);

eig_val2 = eig_val2(1:5); 
eig_vec2 = [];
for i = 1:5
    ind = index(i);
    temp = V(:,index(i));
    vec = [
        (4*temp(1)-1*temp(2))/(3+2*dx*sqrt(16-eig_val2(i))); 
        temp; 
        (4*temp(79)-1*temp(78))/(3+2*dx*sqrt(16-eig_val2(i)))
        %(-4*temp(79)+1*temp(78))/(-3-2*dx*sqrt(16-eig_vals(i)))
    ];
    
    eig_vec2 = [eig_vec2 vec];
end
eig_vec2
eig_vec2 = abs(eig_vec2./sqrt(trapz(xspan, eig_vec2.^2)));
% c = trapz(xspan, vec.^2, 1)
% norm_vec = vec/sqrt(c)

% eig_vals = e(1:5);
% eig_vals
% eig_vec = sort(eig_vec, 2);



% figure 
% hold on
% 
% plot(xspan, eig_vec(:,1))
% plot(xspan, eig_vec(:,2))
% plot(xspan, eig_vec(:,3))
% plot(xspan, eig_vec(:,4))
% plot(xspan, eig_vec(:,5))



A7 = eig_vec2(:,1);
A8 = eig_vec2(:,2);
A9 = eig_vec2(:,3);
A10 = eig_vec2(:,4);
A11 = eig_vec2(:,5);
A12 = eig_val2;

save('A7.dat','A7','-ascii');
save('A8.dat','A8','-ascii');
save('A9.dat','A9','-ascii');
save('A10.dat','A10','-ascii');
save('A11.dat','A11','-ascii');
save('A12.dat','A12','-ascii');

%%

clear all; close all; clc

L = 2;
gamma = 0.05;
xspan = [-L:0.1:L];

A = 0.1;
eps_start = 0; 

k = 1;
tol = 1e-4;

eig_vals = [];
eig_vecs = [];

for modes = 1:2
    eps = eps_start;
    deps = 0.1;
    
    for i = 1:1000
        yinit = [A; A*sqrt(k*L^2-eps)];
        % yinit
        %options = odeset('RelTol',1e12,'AbsTol',1e12);
        [x,y] = ode45(@(x,y) bvpfunc2(x,y,k,gamma,eps), xspan, yinit);
        y1 = y(:,1);
        y2 = y(:,2);
        % calculate norm
        ynorm = trapz(xspan, y(:,1).^2);
        
        
        if abs(ynorm - 1) < tol
            % eig_vals = [eig_vals; eps];
            % eig_vecs = [eig_vecs abs(y(:,1))];
            break;
        else
            A = A/sqrt(ynorm);
        end
        
        % eps shooting
        temp = y(end,2)+sqrt(L^2-eps)*y(end,1);
        
        if abs(temp) < tol 
            % eig_vals = [eig_vals; eps];
            % y_norm = y(:,1)./sqrt(trapz(xspan, y(:,1).^2));
            % eig_vecs = [eig_vecs abs(y_norm)];
            break;
        end
        
        if (-1)^(modes+1)*temp > 0
            eps = eps + deps;
        else
            eps = eps - deps/2;
            deps = deps/2;
        end
    end
    eig_vals = [eig_vals; eps];
    eig_vecs = [eig_vecs y(:,1)];
    eps_start = eps + 0.1;
end

A13 = abs(eig_vecs(:,1));
A14 = abs(eig_vecs(:,2));
A15 = eig_vals;

save('A13.dat','A13','-ascii');
save('A14.dat','A14','-ascii');
save('A15.dat','A15','-ascii');

%% 
% clear variables; close all; clc



clear all; close all; clc

L = 2;
gamma = -0.05;
xspan = [-L:0.1:L];

A = 0.1;
eps_start = 0; 

k = 1;
tol = 1e-4;

eig_vals = [];
eig_vecs = [];

for modes = 1:2
    eps = eps_start;
    deps = 0.1;
    
    for i = 1:1000
        yinit = [A; A*sqrt(k*L^2-eps)];
        % yinit
        %options = odeset('RelTol',1e12,'AbsTol',1e12);
        [x,y] = ode45(@(x,y) bvpfunc2(x,y,k,gamma,eps), xspan, yinit);
        y1 = y(:,1);
        y2 = y(:,2);
        % calculate norm
        ynorm = trapz(xspan, y(:,1).^2);
        
        
        if abs(ynorm - 1) < tol
            % eig_vals = [eig_vals; eps];
            % eig_vecs = [eig_vecs abs(y(:,1))];
            break;
        else
            A = A/sqrt(ynorm);
        end
        
        % eps shooting
        temp = y(end,2)+sqrt(L^2-eps)*y(end,1);
        
        if abs(temp) < tol 
            % eig_vals = [eig_vals; eps];
            % y_norm = y(:,1)./sqrt(trapz(xspan, y(:,1).^2));
            % eig_vecs = [eig_vecs abs(y_norm)];
            break;
        end
        
        if (-1)^(modes+1)*temp > 0
            eps = eps + deps;
        else
            eps = eps - deps/2;
            deps = deps/2;
        end
    end
    eig_vals = [eig_vals; eps];
    eig_vecs = [eig_vecs y(:,1)];
    eps_start = eps + 0.1;
end

A16 = abs(eig_vecs(:,1));
A17 = abs(eig_vecs(:,2));
A18 = eig_vals;

save('A16.dat','A16','-ascii');
save('A17.dat','A17','-ascii');
save('A18.dat','A18','-ascii');