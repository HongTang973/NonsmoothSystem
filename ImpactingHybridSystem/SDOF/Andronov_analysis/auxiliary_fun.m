clc
clear
close all
mu=-0.1;
x_0=-1;
y_0=0.1;
x=@(tau)  exp(-mu*tau)*(x_0*cos(tau)+(y_0*sqrt(1+mu^2)+mu*x_0)*sin(tau));
y=@(tau)  exp(-mu*tau)*(y_0*cos(tau)-(x_0*sqrt(1+mu^2)+mu*y_0)*sin(tau));
F=@(tau) 1-exp(mu*tau)*(cos(tau)-mu*sin(tau));
phi=@(tau,mu) 1-exp(mu*tau)*(cos(tau)-mu*sin(tau));

% set preparation for root finding
    options = optimoptions('fsolve','Display','off','FunctionTolerance',1e-9);

t0=1.5*pi;
[dt1,fval] = fsolve(F,t0,options);
dt1/pi
V=-exp(mu*dt1)*phi(dt1,-mu)/sqrt(1+mu^2)/sin(dt1)
V_=exp(-mu*dt1)*phi(dt1,mu)/sqrt(1+mu^2)/sin(dt1)
tau_list=-2*pi:0.01:2*pi;

% phi
phi_list=zeros(length(tau_list),1);
Dis=zeros(length(tau_list),1);
Vt=zeros(length(tau_list),1);
for i=1:length(tau_list)
    phi_list(i)=phi(tau_list(i),mu); 
%     phi_list(i)=exp(-2*mu*tau_list(i))*phi(tau_list(i),mu)/phi(tau_list(i),-mu); 
     Dis(i)=x(tau_list(i)); 
      Vt(i)=y(tau_list(i)); 
end
figure
plot(tau_list,phi_list,'k-','linewidth',1.2)
set(gca,'fontname','times new roman','fontsize',12)
figure
plot(tau_list,Dis,'k-','linewidth',1.2)
set(gca,'fontname','times new roman','fontsize',12)
figure
plot(tau_list,Vt,'k-','linewidth',1.2)
set(gca,'fontname','times new roman','fontsize',12)


tau_list2=2*pi:-0.01:pi;
dSS=zeros(length(tau_list2),1);
for i=1:length(tau_list2)
    dSS(i)=phi(tau_list2(i),-mu)/phi(tau_list2(i),mu); 
end
figure
plot(tau_list2,dSS,'k-','linewidth',1.2)
set(gca,'fontname','times new roman','fontsize',12)