clc
clear
close all
mu=20;
x_0=-1;
y_0=0.1;
x=@(tau)  exp(-mu*tau).*(x_0*cos(tau)+(y_0*sqrt(1+mu^2)+mu*x_0).*sin(tau));
y=@(tau)  exp(-mu*tau).*(y_0*cos(tau)-(x_0*sqrt(1+mu^2)+mu*y_0).*sin(tau));
F=@(tau) 1-exp(mu*tau).*(cos(tau)-mu*sin(tau));
%> define the auxilary function
phi=@(tau,mu) 1-exp(mu*tau).*(cos(tau)-mu*sin(tau));
V_out = @(tau,mu) exp(-mu*tau).*phi(tau,mu)./sin(tau);

% set preparation for root finding
options = optimoptions('fsolve','Display','off','FunctionTolerance',1e-9);
t0=1.5*pi;
[dt1,fval] = fsolve(F,t0,options);
% dt1/pi
V=-exp(mu*dt1)*phi(dt1,-mu)/sqrt(1+mu^2)/sin(dt1)
V_=exp(-mu*dt1)*phi(dt1,mu)/sqrt(1+mu^2)/sin(dt1)

tau_list=0:0.01:2*pi;

% phi
% phi_list=zeros(length(tau_list),1);
% Dis=zeros(length(tau_list),1);
% Vt=zeros(length(tau_list),1);
phi_list=phi(tau_list,mu);
phi_m_list = phi(tau_list,-mu);

V_list  = V_out(tau_list,mu);
Dis=x(tau_list);
Vt=y(tau_list);

figure
plot(tau_list,phi_list,'k-','linewidth',1.2)
hold on
plot(tau_list,phi_m_list,'k-','linewidth',1.2)
plot(dt1,fval,'ro')
plot([-2*pi 2*pi],[0 0],'k-')
set(gca,'fontname','times new roman','fontsize',12)
ylim([-1 5])
figure
plot(tau_list,V_list)
hold on
plot([-2*pi 2*pi],[0 0],'k-')
figure
plot(tau_list,Dis,'k-','linewidth',1.2)
set(gca,'fontname','times new roman','fontsize',12)
figure
plot(tau_list,Vt,'k-','linewidth',1.2)
set(gca,'fontname','times new roman','fontsize',12)


% tau_list2=2*pi:-0.01:pi;
% dSS=zeros(length(tau_list2),1);
% for i=1:length(tau_list2)
%     dSS(i)=phi(tau_list2(i),-mu)/phi(tau_list2(i),mu); 
% end
% figure
% plot(tau_list2,dSS,'k-','linewidth',1.2)
% set(gca,'fontname','times new roman','fontsize',12)