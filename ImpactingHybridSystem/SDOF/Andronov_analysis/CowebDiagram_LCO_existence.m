clc
clear
close all
% SDOF the LCO existence condition
mu=0.2;

phi=@(tau,mu) 1-exp(mu*tau).*(cos(tau)-mu*sin(tau));

% set preparation for root finding
F=@(tau) 1-exp(mu*tau)*(cos(tau)-mu*sin(tau));
options = optimoptions('fsolve','Display','off','FunctionTolerance',1e-9);
% get the critical V
t0=1.5*pi;
[dt1,fval] = fsolve(F,t0,options);
dt1/pi;
% dt1=pi;

%  V=-exp(mu*dt1).*phi(dt1,-mu)/(sqrt(1+mu^2))./sin(dt1)

tau_list=dt1:-0.00001:pi+0.01;
V=-exp(mu*tau_list).*phi(tau_list,-mu)./(sqrt(1+mu^2))./sin(tau_list);
FV=exp(-mu*tau_list).*phi(tau_list,mu)./(sqrt(1+mu^2))./sin(tau_list);
%
r=2;
V_=0:0.1:max(V);
PV=V_./r;
figure
plot(V/100,-FV)
hold on
plot(V_/100,PV)

cline=log(abs(FV)./V)/pi;


