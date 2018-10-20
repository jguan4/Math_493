clear all
%close all

p=2; %number of parameters
nsteps = 5000;
burntime = 100;
D=0.01;

load('.\data4\data4.mat');
t=xx;
sz = length(t);
sigma=0.03;
syn_data = yy;

%Spring model
solve=@(C,beta)(C*exp(-beta*t));

theta1 = zeros(1,p);
%  initial guess for p parameters


y=solve(theta1(1),theta1(2));
% solve  gives the vector y at the desired time points t for parameter theta1


chi1= sum((syn_data - y).^2);

for n = 1: nsteps
    
    theta2 =theta1 + D.* randn(1,p); %this calculates proposed new parameter value
    % p is the number of parameters and
    % D is the guess standard deviation
    
    y=solve(theta2(1),theta2(2)); % model response for the new parameter value
    
    chi2=sum((syn_data-y).^2); %SS error for the new value
    ratio=exp((-chi2+chi1)/(2*sigma^2)) %ratio r
    
    if rand < ratio % if ratio>1, holds automatically, otherwise holds with probability r=ratio
        theta1=theta2;
        chi1=chi2;
    end
    
    if n > burntime
        thetasave(n-burntime,:)=theta1;
        chisave(n-burntime,:)=chi1;
    end
    
end

figure(1);

subplot(1,2,1);
hold on;
axis square;
set(gca,'fontsize', 15);
xlim([1 n-burntime]);
title('Parameter');
plot(thetasave);

subplot(1,2,2);
hold on;
axis square;
set(gca,'fontsize', 15);
xlim([1 n-burntime]);
title('SSE');
plot(chisave);

figure;
hold on;
axis square;
set(gca,'fontsize', 15);
title('Data vs. Estimate')
plot(xx,solve(thetasave(end,1),thetasave(end,2)))
plot(xx,yy)
legend('Estimate','Data')
xlabel('t')
ylabel('x')
drawnow;





