% Test example for nonlinear least squares fitting
% Expected solution: a = 2;

close all;
clear all;
plotting=1;

global xexp
global texp

%Experimental data
xexp = [1.1,-0.98,-0.3]';
texp = [0:2];

method =2; % 1=fminsearch, 2=Newton

if (method==1)
    a_opt = 1;
    % Finding the parameter values minimizing SSE error using built-in function
    [a_opt,fval,exitflag,output] = fminsearch(@SSE,0)
    
else
    %Newton-Raphson method
    aa(1,:)=[0.5 0];
    delta1 = 1e-10;
    delta2 = 1e-10;
    Nsteps = 50;
    
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
        g = J1'*(xexp(1,:)' - z(:,1));
        
        
        %Newton direction = inv(H)*g, but
        %Better to use slash command for a more stable calculation
        aa(k+1,:) = aa(k,:) + (H\g);
        a_opt = aa(k+1,:);
        
        if (norm(aa(k+1)-aa(k))<delta1) | (abs(SSE(aa(k+1))-SSE(aa(k)))<delta2)
            break;
        end
    end
    
    figure(10);
    plot(1:length(aa), aa,'ro');
    
end


if (plotting==1)
    
    xt0 =a_opt(2);
    [tspan,x] = ode45(@(t,x)system(t,x,a_opt),texp,xt0);
    
    figure(1);
    set(gca,'FontName','Arial','FontSize',14,'FontWeight','Bold','LineWidth', 1);
    hold on;
    axis square;
    title('y(1) vs. t');
    plot(tspan,x,'r.');
    plot(texp,xexp,'bo');
    legend('ODE','data');
    
end

function z = SSE(a)
global texp
global xexp

xt0 =[1];
[tspan,x] = ode45(@(t,x)system(t,x,a),texp,xt0);

z = norm(x-xexp,2);

end

function dxdt=system(t,x,a)

dxdt = a*x^2;

end

function dzdt = sens(t,z,a)

dzdt(1)=a*z(1)^2;
dzdt(2)=z(1)^2+2*a*z(1)*z(2);
dzdt(3)=2*a*z(1)*z(3)
dzdt = dzdt';
end
