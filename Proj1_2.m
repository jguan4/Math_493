close all;
clear all;
format long;
plotting=1;

global xexp
global texp

%Experimental data
xexp = [-0.9 -0.3 -0.28 -0.15 -0.13]';
texp = [1:5];
a_opt_init = 1;
xt0_init = -1;

method =2; % 1=fminsearch, 2=Newton

if (method==1)
    % Finding the parameter values minimizing SSE error using built-in function
    [a_opt,fval,exitflag,output] = fminsearch(@SSE,xt0_init)
    
else
    %Newton-Raphson method
    aa(1,:)=[a_opt_init xt0_init];
    error(1) = SSE(aa(1,:));
    delta1 = 1e-16;
    delta2 = 1e-16;
    Nsteps = 200;
    
    for k=1:Nsteps
        
        %Need to solve ODE on a finer grid
        tt = linspace(texp(1),texp(end),100);
        zt0 =[aa(k,2);0;1]; %initial data for the sensitivity system
        [tspan,zz] = ode45(@(t,z)sens(t,z,aa(k,1)),tt,zt0);
        z = interp1(tt,zz,texp); %getting back to the coarser grid to compare with xexp
        
        %Jacobian calculation
        J1 = [z(:,2:3)];
        
        %Hessian with weights
        H = J1'*J1;
        g = J1'*(xexp - z(:,1));
        
        
        %Newton direction = inv(H)*g, but
        %Better to use slash command for a more stable calculation
        aa(k+1,:) = aa(k,:) + (H\g)';
        a_opt = aa(k+1,:);
        error(k+1) = SSE(a_opt);
        
        if (norm(aa(k+1,:)-aa(k,:))<delta1) | (abs(SSE(aa(k+1,:))-SSE(aa(k,:)))<delta2)
            break;
        end
    end
end

% Plotting convergence of parameters
figure(2);
subplot(1,3,1);
hold on;
title('\theta')
axis square;
plot(1:length(aa(:,1)),aa(:,1),'r.-');

subplot(1,3,2);
hold on;
title('x_0')
axis square;
plot(1:length(aa(:,2)),aa(:,2),'r.-');

subplot(1,3,3);
hold on;
title('SSE');
axis square;
plot(1:length(aa(:,1)),error,'k.-');

if (plotting==1)
    
    [tspan,x] = ode45(@(t,x)system(t,x,a_opt(1)),texp,a_opt(2));
    
    figure(1);
    set(gca,'FontName','Arial','FontSize',14,'FontWeight','Bold','LineWidth', 1);
    hold on;
    axis square;
    title('x vs. t');
    plot(tspan,x,'r.');
    plot(texp,xexp,'bo');
    legend('ODE','data');
    
end

function z = SSE(a)
global texp
global xexp

xt0 =a(2);
[tspan,x] = ode45(@(t,x)system(t,x,a(1)),texp,xt0);

z = norm(x-xexp,2);

end

function dxdt=system(t,x,a)

dxdt = a*x^2;

end

function dzdt = sens(t,z,a)

dzdt(1)=a*z(1)^2;
dzdt(2)=z(1)^2+2*a*z(1)*z(2);
dzdt(3)=2*a*z(1)*z(3);
dzdt = dzdt';
end
