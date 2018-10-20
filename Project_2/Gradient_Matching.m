close all;
clear all;
format long;

load('.\data4\data4.mat');
t=xx;
syn_data = yy;

num_poly = 11; % number of polynomial approximation in collocation
k = 5; % order of collocation
n_steps = 100;
tol = 1e-3;
theta_est(1)=1; % parameter estimation initialization

tau = linspace(t(1),t(end),201);
knots = augknt(linspace(t(1),t(end),num_poly),k);
colmat = spcol(knots,k,brk2knt(tau,2));
phi = colmat(1:2:end,:);
phi_t = phi(1:5:end,:);
Dphi =colmat(2:2:end,:);
Dphi_t = Dphi(1:5:end,:);
coefs = ((phi_t')*phi_t)\(phi_t'*syn_data');
x_est = phi_t*coefs;
xt_0 = x_est(1);
Dx_est = Dphi_t*coefs;

for n = 2: n_steps
    J = x_est;
    H = x_est'*x_est;
    g = x_est'*(Dx_est+theta_est(n-1).*x_est);
    theta_est(n) = theta_est(n-1)-H\g;
    
    error(n-1) = SSE(t,xt_0, theta_est(n), syn_data);
    if (theta_est(n)-theta_est(n-1))^2<1e-16 && error(n-1)<tol
        break;
    end
end

figure;
hold on;
axis square;
set(gca,'fontsize', 15);
xlim([1 n]);
title('Parameter');
plot(theta_est);

figure
hold on;
axis square;
set(gca,'fontsize', 15);
xlim([1 n]);
title('SSE');
plot(error);

[tspan,x] = ode45(@(t,x)system(t,x,theta_est(end)),t,xt_0);

figure;
set(gca,'FontName','Arial','FontSize',14,'FontWeight','Bold','LineWidth', 1);
hold on;
axis square;
title('x vs. t');
plot(tspan,x,'r');
plot(t,syn_data,'b');
legend('ODE','data');
drawnow;

function z = SSE(t,xt_0, theta, syn_data)

[tspan,x] = ode45(@(t,x)system(t,x,theta),t,xt_0);

z = norm(x-syn_data,2);

end

function dxdt=system(t,x,a)

dxdt = -a*x;

end