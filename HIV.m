format long
global tspan
global texp
global xexp
global xt0

close all;

load("/home/meglstudent/Downloads/hiv_data.mat");
tspan = (0:.1:200);
texp = hiv_data(:,1);
xexp = hiv_data(:,2:end);
xt0 = [.9e6 4000 .1 .1 1 12];

mcmc_flag=0; % 1 to do mcmc, 0 to skip and just graph

n=length(xexp); %# of data points
p=6; % #of parameters
nsteps = 10000;
D=.03;
burntime=1;

q0 = [.3 .7 .01 1e-4 1e4 100]; %initial parameter guesses can be based on fminsearch etc

sigmas = [750, 20, 100, 15, 1000, 1]; %These are guesses at the variance of the variables in an effort to limit sse size 

if exist('Q')==0
    Q = fminsearch(@(Q)SSE (Q,tspan,xt0,xexp,texp),q0);
end

 %system solution using the found parameters
[t,y] = ode45(@(t,X) sys(t,X,Q),tspan,xt0);

theta1 = Q;
%  initial guess for p parameters

% solve  gives the vector y at the desired time points t for parameter theta1

y=interp1(t,y,texp);

chi1 = sum(sum((xexp - y).^2).*(1./(2*sigmas.^2)));

if mcmc_flag == 1
    thetasave = zeros(nsteps,p);
    chisave = zeros(nsteps);
    for n = 1: nsteps+1
        if mod(n,25)==0
            n
        end
        theta2 = abs(theta1 + randn(1,p).*(Q./25)); %this calculates proposed new parameter value
        %theta2
        %theta1
        % p is the number of parameters and
        % D is the guess standard deviation

        [t,y] = ode45(@(t,X) sys(t,X,theta2),tspan,xt0); % model response for the new parameter value
        y=interp1(t,y,texp);

        chi2=sum(sum((xexp - y).^2).*(1./(sigmas.^2))); %SS error for the new value
        %chi2
        %chi1
        ratio=exp((-chi2+chi1));%/(2*sigma^2)%); %ratio r
        %ratio

        if rand < ratio % if ratio>1, holds automatically, otherwise holds with probability r=ratio
            theta1=theta2;
            chi1=chi2;
        end

        if n > burntime
            thetasave(n-burntime,:)=theta2;
            chisave(n-burntime)=chi1;
        end

    end
end

means = mean(thetasave);
[t,y] = ode45(@(t,X) sys(t,X,means),tspan,xt0);
y=interp1(t,y,texp);


stddevs = std(y);

xhat= mean(texp);
Sxx=0;
for i=[1:length(texp)]
    Sxx=Sxx+(texp(i)-xhat)^2;
end
upsilon = zeros(length(xexp),p);
for i=[1:length(texp)]
    for j = [1:p]
        upsilon(i,j) = 1.96 * stddevs(j) *(1/length(texp) + (texp(i)-xhat)^2/Sxx)^(1/2);
    end
end

ub = y+upsilon;
lb = y-upsilon;
figure(4);
for i = [1:length(xt0)]
    subplot(3,2,i)
    hold on;
    plot(texp,ub(:,i));
    plot(texp,lb(:,i));
    plot(texp,y(:,i));
    plot(texp,xexp(:,i));
end



figure(1);

subplot(1,2,1);
hold on;
axis square;
set(gca,'fontsize', 15);
%axis([1 n-burntime 0.5 2.5]);
title('MCMC-estimated parameter');
plot(thetasave(:,1));

subplot(1,2,2);
hold on;
axis square;
set(gca,'fontsize', 15);
%axis([1 n-burntime 0 0.03]);
title('Sum of squares error');
plot(chisave);
drawnow;

%density distributions
figure(2);
var_names = {'$b_{E}$' '$\delta$' '$d_{1}$' '$k_{2}$' '$\lambda_{1}$' '$K_{b}$'};
for i = [1:6]
    subplot(3,2,i);
    hold on;
    title(var_names(i),'interpreter','latex');
    h = histc(thetasave(burntime:end,i),(min(thetasave(burntime:end,i)):(max(thetasave(burntime:end,i))-min(thetasave(burntime:end,i)))/25:max(thetasave(burntime:end,i))));
    plot(h);
    xticklabels(min(thetasave(burntime:end,i)):(max(thetasave(burntime:end,i))-min(thetasave(burntime:end,i)))/25:max(thetasave(burntime:end,i)));
   
    
end

%Correlation Scatter Plots
figure(3);
for i=[1:5]
    for j=[2:6]
        if i<j
            hold on;
            subplot(5,5,(i-1)*5+(j-1));
            scatter(thetasave(burntime:end,i),thetasave(burntime:end,j));
            lsline;
            xlabel(var_names(i),'interpreter','latex');
            ylabel(var_names(j),'interpreter','latex');
        end
    end
end




function sse = SSE(Q,tspan,xt0,xexp,texp)
    for i = [1:length(Q)]
        if Q(i)<0
            sse=inf;
        else
            [t,x] = ode45(@(t,X) sys(t,X,Q),tspan,xt0);
            x = interp1(t,x,texp);
            sse=sum(sum((xexp-x)'*(xexp-x)));
        end
    end
end

function dxdt = sys(t,X,Q)
    be=Q(1); delta=Q(2); d1=Q(3); k2=Q(4); lambda1=Q(5); Kb=Q(6);
    d2 = 0.01;
    k1 = 8e-7;
    lambda2 = 31.98;
    eps = .1;
    de = 0.25;
    Nt = 100;
    f = 0.34;
    m1 = 1e-5;
    m2 = 1e-5;
    lambda_e = 1;
    delta_e = 0.1;
    c = 13;
    p1 = 1;
    p2 = 1;
    Kd = 1;
    T1 = X(1);
    T2 = X(2);
    T1s = X(3);
    T2s = X(4);
    V = X(5);
    E = X(6);
   
    dT1 = lambda1 - d1*T1 - (1-eps)*k1*V*T1;
    dT2 = lambda2 - d2*T2 - (1-f*eps)*k2*V*T2;
    dT1s = (1-eps)*k1*V*T1 - delta*T1s - m1*E*T1s;
    dT2s = (1 - f*eps)*k2*V*T2 - delta*T2s - m2*E*T2s;
    dV = Nt*delta*(T1s + T2s) - c*V -((1-eps)*p1*k1*T1 + (1-f*eps)*p2*k2*T2)*V;
    dE = lambda_e + E*(be*(T1s + T2s))/(T1s + T2s + Kb) - E*(de*(T1s + T2s))/(T1s + T2s + Kd) - delta_e*E;
   
    dxdt = [dT1;dT2;dT1s;dT2s;dV;dE];
end
