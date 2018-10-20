close all;
clear all;
format long;
% [3.12912405083164,5.34662575325762]

load('.\data4\data4.mat');
t=xx;
N=length(t);
dt = t(2)-t(1);
syn_data = yy;

num_poly = 11; % number of polynomial approximation in collocation
k = 5; % order of collocation
pen_p=0.5;
n_steps = 100;
tol = 1e-10;
theta_est(1)=1; % parameter estimation initialization
coefs_est(1,:) = ones([1,num_poly+k-2]);

tau = linspace(t(1),t(end),201); 
knots = augknt(linspace(t(1),t(end),num_poly),k);
colmat = spcol(knots,k,brk2knt(tau,2));

for n = 2: n_steps
    % calculate coef
    phi = colmat(1:2:end,:);
    phi_t = phi(1:5:end,:);
    Dphi =colmat(2:2:end,:);
    Dphi_t = Dphi(1:5:end,:);
    R_theta = calculateR(Dphi,phi,theta_est(n-1),tau(2)-tau(1));
    A = (1-pen_p)*(phi_t')*phi_t/N + pen_p*R_theta/t(end);
    coefs_est(n,:) = (A\((1-pen_p).*(phi_t'*syn_data')/N))';
    new_approx_data = phi_t*coefs_est(n,:)';
    
    % calculate theta
    dAdtheta = calculatedAdtheta(pen_p, t(end), tau(2)-tau(1), Dphi, phi, theta_est(n-1));
    A_inv = inv(A);
    dMdtheta = (1-pen_p)*(phi_t)*(-A_inv*dAdtheta*A_inv)*(phi_t')/N;
    dx_est_dtheta = dMdtheta * syn_data';
    H = dx_est_dtheta'*dx_est_dtheta;
    g = dx_est_dtheta'*(syn_data'-new_approx_data);
    theta_est(n) = theta_est(n-1)+H\g;
    
    error(n-1) = (theta_est(n)-theta_est(n-1))^2 + norm((coefs_est(n,:)-coefs_est(n-1,:)),2);
    if error(n-1) < tol
        break;
    end
end

figure;
subplot(1,2,1);
hold on;
axis square;
set(gca,'fontsize', 15);
xlim([1 n]);
title('Parameter');
plot(theta_est);

subplot(1,2,2);
hold on;
axis square;
set(gca,'fontsize', 15);
xlim([1 n]);
title('SSE');
plot(error);

figure;
hold on;
axis square;
set(gca,'fontsize', 15);
title('Data vs. Estimate')
plot(t,phi_t*coefs_est(end,:)')
plot(t,syn_data)
legend('Estimate','Data')
xlabel('t')
ylabel('x')
drawnow;


function R = calculateR(Dphi, phi, beta,dt)
R = zeros([size(phi,2),size(phi,2)]);

for k=1:size(phi,2)
    for l=1:size(phi,2)
        for tt = 1:size(phi,1)-1
            R(k,l) = R(k,l) + dt*(Dphi(tt,k) + beta.* phi(tt,k))*(Dphi(tt,l) + beta.* phi(tt,l));
        end
    end
end

end

function dAdtheta = calculatedAdtheta(pen_p, T, dx, Dphi, phi, beta)
dAdtheta = zeros([size(phi,2),size(phi,2)]);

for k=1:size(phi,2)
    for l=1:size(phi,2)
        for tt = 1:size(phi,1)-1
            dAdtheta(k,l) = dAdtheta(k,l) + (Dphi(tt,k) + beta* phi(tt,k))*(phi(tt,l)) + ...
                (Dphi(tt,l) + beta* phi(tt,l))*(phi(tt,k));
        end
    end
end

dAdtheta = pen_p*dx *dAdtheta/T;
end