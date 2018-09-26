close all;
clear all;
format long;
plotting=1;

global xexp
global texp
global xt0

%Experimental data
xexp = [0.9 0.858 0.78 0.7 0.6 0.5;
    0.1 0.14 0.2 0.25 0.3 0.38]';
texp = [0:5];
xt0 = [0.9 0.1];

method =2; % 1=fminsearch, 2=Newton

if (method==1)
    % Finding the parameter values minimizing SSE error using built-in function
    [a_opt,fval,exitflag,output] = fminsearch(@SSE,[0 0])
    
else
    %Newton-Raphson method
    aa(1,:)=[0 0];
    error(1) = SSE(aa(1,:));
    delta1 = 1e-16;
    delta2 = 1e-16;
    Nsteps = 200;
    
    for k=1:Nsteps
        
        %Need to solve ODE on a finer grid
        tt = linspace(texp(1),texp(end),300);
        zt0 =[xt0 0 0 0 0]; %initial data for the sensitivity system
        [tspan,zz] = ode45(@(t,z)sens(t,z,aa(k,:)),tt,zt0);
        z = interp1(tt,zz,texp); %getting back to the coarser grid to compare with xexp
        
        %Jacobian calculation
        J1 = [z(:,3:4)];
        J2 = [z(:,5:6)];
        
        %Hessian with weights
        H = J1'*J1 + J2'*J2;
        g = J1'*(xexp(:,1) - z(:,1)) + J2'*(xexp(:,2) - z(:,2));
        
        
        %Newton direction = inv(H)*g, but
        %Better to use slash command for a more stable calculation
        aa(k+1,:) = aa(k,:) + (H\g)';
        a_opt = aa(k+1,:);
        error(k+1) = SSE(a_opt);
        
        if (norm(aa(k+1,:)-aa(k,:))<delta1) | (abs(SSE(aa(k+1,:))-SSE(aa(k,:)))<delta2)
            tt=[0:100];
            [tspan,xd] = ode45(@(t,x)system(t,x,a_opt),tt,xt0);
            figure(3)
            set(gca,'FontName','Arial','FontSize',14,'FontWeight','Bold','LineWidth', 1);
            hold on;
            xlabel("Time (days)")
            ylabel("Population Percentage")
            axis square;
            plot(tspan,xd,'LineWidth',2)
            legend("Susceptible","Infected")
            
            % Plotting convergence of parameters
            figure(2);
            subplot(1,3,1);
            set(gca,'FontName','Arial','FontSize',14,'FontWeight','Bold','LineWidth', 1);
            hold on;
            title('a')
            xlabel("Iteration Step")
            axis square;
            plot(1:length(aa(:,1)),aa(:,1),'r.-');
            
            subplot(1,3,2);
            set(gca,'FontName','Arial','FontSize',14,'FontWeight','Bold','LineWidth', 1);
            hold on;
            title('b')
            xlabel("Iteration Step")
            axis square;
            plot(1:length(aa(:,2)),aa(:,2),'r.-');
            
            subplot(1,3,3);
            set(gca,'FontName','Arial','FontSize',14,'FontWeight','Bold','LineWidth', 1);
            hold on;
            title('SSE');
            xlabel("Iteration Step")
            axis square;
            plot(1:length(aa(:,1)),error,'k.-');
            break;
        end
    end
end

if (plotting==1)
    
    [tspan,x] = ode45(@(t,x)system(t,x,a_opt),texp,xt0);
    
    figure(1);
    subplot(1,2,1);
    set(gca,'FontName','Arial','FontSize',14,'FontWeight','Bold','LineWidth', 1);
    hold on;
    axis square;
    title('x vs. t');
    plot(tspan(:,1),x(:,1),'r.');
    plot(texp,xexp(:,1),'bo');
    legend('ODE','data');
    
    subplot(1,2,2);
    set(gca,'FontName','Arial','FontSize',14,'FontWeight','Bold','LineWidth', 1);
    hold on;
    axis square;
    title('y vs. t');
    plot(tspan,x(:,2),'r.');
    plot(texp,xexp(:,2),'bo');
    legend('ODE','data');
    
end

function z = SSE(a)
global texp
global xexp
global xt0

[tspan,x] = ode45(@(t,x)system(t,x,a),texp,xt0);

z = norm(x-xexp,2);

end

function dxdt=system(t,x,a)

dxdt(1) = -a(1)*x(1)*x(2);
dxdt(2) = a(1)*x(1)*x(2) - a(2)*x(2);

dxdt = dxdt';

end


function dzdt = sens(t,z,a)
dzdt(1) = -a(1)*z(1)*z(2);
dzdt(2) = a(1)*z(1)*z(2)-a(2)*z(2);
dzdt(3) = -z(1)*z(2) - a(1)*z(2)*z(3) - a(1)*z(2)*z(4);
dzdt(4) = z(1)*z(2) + a(1)*z(2)*z(3) + (a(1)*z(1)-a(2))*z(4);
dzdt(5) = -a(1)*z(2)*z(5) -a(1)*z(1)*z(6);
dzdt(6) = -z(2) + a(1)*z(2)*z(5) + (a(1)*z(1)-a(2))*z(6);
dzdt = dzdt';
end