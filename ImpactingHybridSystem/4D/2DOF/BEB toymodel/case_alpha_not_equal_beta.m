% alpha not equal lambda
clc
close all
theta=pi/6;
alpha = -0.1;
beta  = 0.2;
lambda = -0.2;
gamma =(alpha-lambda)/beta;
phi_1 = @(tau) 1-exp(gamma*tau).*cos(tau);
phi_2 = @(tau) phi_1(tau) +gamma*exp(gamma*tau).*sin(tau);
phi_3 = @(tau) phi_1(tau) - exp(tau*gamma).*sin(tau)/gamma;


tau=0:0.01:2*pi;
p1 = phi_1(tau);
p2 = phi_2(tau);
p3 = phi_3(tau);
p0 = 0*tau;


figure(1)
plot(tau/pi,p1,'r-','displayname','$\varphi_1(\tau)$')
hold on
plot(tau/pi,p0,'b-','displayname','$\rm{zero~line}$')
plot(tau/pi,p2,'k-','displayname','$\varphi_2(\tau)$')
% plot(tau/pi,p3,'g-','displayname','$\varphi_3(\tau)$')
xlabel('$\tau(\pi)$','interpreter','latex')
legend('interpreter','latex')
% detect the region of the y0 z0
% obeserve the z0 and y0

C=1;
AA=C-lambda;
y0=(-gamma*beta + gamma*beta*exp(-lambda*tau/beta) - AA*exp(tau*gamma).*cos(tau) + AA)./(phi_2(tau)*cos(theta)*beta);
z0 =cot(theta) + (beta*(exp(-lambda*tau/beta) - 1) - AA*exp(tau*gamma).*sin(tau))./(beta*sin(theta)*cos(theta).*phi_2(tau));

figure
plot(tau/pi,y0,'r-','displayname','$y_0$')
hold on
plot(tau/pi,z0,'k-','displayname','$z_0$')
xlabel('$\tau(\pi)$','interpreter','latex')
legend('interpreter','latex')

%
% ratio of the poincare
r=(AA*exp(tau*(2*beta*gamma + lambda)/beta)+...
    ((beta*gamma^2 - (AA + lambda)*gamma + lambda*gamma + beta)*sin(tau) - cos(tau)*AA).*exp((beta*gamma + lambda)*tau/beta)+ ...
    ((-beta*gamma^2 - gamma*lambda - beta)*sin(tau) + lambda*cos(tau)).*exp(tau*gamma) - lambda)./(phi_2(tau)*(AA + lambda));

figure
plot(tau/pi,r,'s-','displayname','$ratio_r(\tau)$')
legend('interpreter','latex')
% define ratio_Z
D1 = AA*exp(lambda*tau/beta).*(cos(theta)^2*gamma*(-phi_1(tau)) + gamma*phi_3(tau));
D2 = beta*(gamma^2 + 1)*exp(lambda*tau/beta).*(exp(-lambda*tau/beta) - sin(theta)^2);
D3 = -beta*cos(theta)^2*(exp(lambda*tau/beta).*(1 - phi_2(tau)) + gamma^2);
N1 = -AA*exp(2*(gamma + lambda/beta)*tau).*(cos(theta)^2*gamma*(-phi_1(-tau)) + gamma*phi_3(-tau));
N2 = beta*(gamma^2 + 1)*cos(tau).*exp(lambda*tau/beta).*exp(tau*gamma).*(sin(theta)^2 - exp(lambda*tau/beta));
N3 = beta*exp(lambda*tau/beta)*cos(theta)^2.*(gamma^2*exp(lambda*tau/beta).*(1-phi_3(tau)) + 1);

rz = (N1+N2+N3)./(D1+D2+D3);
Asyp1 = (cos(theta)^2*gamma*(-phi_1(tau)) + gamma*phi_3(tau));
Asyp2 = (cos(theta)^2*gamma*(-phi_1(-tau)) + gamma*phi_3(-tau));
figure
plot(tau/pi,rz,'p-','displayname','$ratio_Z(\tau)$')
legend('interpreter','latex')
%
figure
plot(tau/pi,r,'s-','displayname','$ratio_r(\tau)$')
legend('interpreter','latex')
xlim([0 1.5])
ylim([0 2])
%
figure
plot(tau/pi,rz,'p-','displayname','$ratio_Z(\tau)$')
legend('interpreter','latex')
xlim([0 1])
ylim([0 2])

% detect the pole of the denom of r_Z
nr_Z = N1+N2+N3;
% dr_Z = D1+D2+D3;
dr_Z = (-cos(tau)*cos(theta)^2*beta - sin(tau)*(C - lambda)).*exp(lambda*tau/beta) + (cos(theta) - 1)*beta*(cos(theta) + 1).*exp(lambda*tau/beta) + beta;
% dr_Z = (((AA*gamma - beta)*cos(tau) + beta*gamma*sin(tau))*cos(theta)^2 ...
%     - (gamma*cos(tau) + sin(tau))*AA).*exp((beta*gamma + lambda)*tau/beta) ...
%     - (-beta*gamma^2 + AA*gamma - beta)*(cos(theta) + 1)*(cos(theta) - 1)*exp(lambda*tau/beta) ...
%     - beta*(cos(theta)^2*gamma^2 - gamma^2 - 1);
figure
plot(tau/pi,dr_Z,'p-','displayname','$dr_ratio_Z(\tau)$')
hold on
plot(tau/pi,nr_Z,'s-','displayname','$nr_ratio_Z(\tau)$')
legend('interpreter','latex')

% 
% gamma= 0:0.01:10;
% half_pi = (gamma.*beta*cos(theta)^2 - AA).*exp(((beta*gamma + lambda)*pi)/(2*beta)) - ...
%     (-beta*gamma.^2 + AA*gamma - beta)*(cos(theta) + 1)*(cos(theta) - 1).*exp(lambda*pi/(2*beta))...
%     - beta*(cos(theta)^2*gamma.^2 - gamma.^2 - 1);
% figure
% plot(gamma,half_pi)

%
dr_C = -exp(tau*(2*beta*gamma + lambda)/beta)*lambda + ((beta*gamma^2 + gamma*lambda + beta)*sin(tau) + cos(tau)*lambda).*exp((beta*gamma + lambda)*tau/beta)...
    + ((-beta*gamma^2 - gamma*lambda - beta)*sin(tau) + cos(tau)*lambda).*exp(tau*gamma) - lambda;
figure
plot(tau/pi,dr_C)
ylim([-2 2])
%
phi_4 =@(tau)((-beta*gamma^2 - gamma*lambda - beta).*sin(tau) + cos(tau)*lambda).*exp(tau*gamma) - lambda;
pminus = phi_4(-tau);
pposit = phi_4(tau);
figure
plot(tau/pi,pminus)
hold on
plot(tau/pi,pposit)
