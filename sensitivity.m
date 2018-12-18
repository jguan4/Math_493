format longg
global tspan
global texp
global xexp
global xt0

load('hiv_data.mat')
tspan = (0:.1:200); %t for solving
texp = hiv_data(:,1); %t data
xexp = hiv_data(:,2:end); %x data
xt0 = [.9e6 4000 .1 .1 1 12]; %initial conditions
 
n=length(xexp); %# of data points
p=6; % #of parameters
 
q0 = [.3 .7 .01 1e-4 1e4 100]; %initial parameter guesses can be based on fminsearch etc

 %system solution using the found parameters
[t,x] = ode45(@(t,X) sys(t,X,q0),tspan,xt0);
xslv = interp1(t,x,texp);

for i = 1:6
        subplot(3,2,i)
        plot(texp,xexp(:,i))
        hold on
        plot(t,x(:,i))
        xt = sprintf('Plot of Equation %d', i);
        title(xt,'interpreter','latex');
        grid
end

%Now we have our mean value for our parameters, we will 

sse = @(dta, slv) diag((dta - slv)' * (dta - slv)); %Returns sums of squared error vector for each parameter

s2 = sse(xexp, xslv)/(n-p); %initial variance estimate for each parameter

%next we will find the sensitivty to each parameter 

%Let sens be the sensitivty tensor

epsilon = .000005;
sens = zeros(n,p,6);
for i = 1:6
    ei = zeros(1,6);
    ei(i) = 1;
    qeps = q0 + epsilon*ei;
    
    [t,soltemp] = ode45(@(t,X) sys(t,X,qeps),tspan,xt0);
    soltemp = interp1(t,soltemp,texp);
    sensitivity_i = (soltemp - xslv)/epsilon;
    sens(:,:,i) = sensitivity_i;
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



    
    

