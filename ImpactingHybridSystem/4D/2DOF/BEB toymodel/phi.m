% check the roots
clc
close all
theta=pi/6;
alpha = -0.5;
beta  = 0.2;
lambda = -0.1;
gamma =(alpha-lambda)/beta;
phi_1 = @(tau) 1-exp(gamma*tau).*cos(tau);
phi_2 = @(tau) phi_1(tau) +gamma*exp(gamma*tau).*sin(tau);
phi_3 = @(tau) phi_1(tau) - exp(tau*gamma).*sin(tau)/gamma;


tau=0:0.01:4*pi;

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

% obeserve the z0 and y0

C=0.2428;
AA=C-lambda;
y0=(beta*gamma*exp(-lambda*tau/beta) - beta*gamma - exp(tau*gamma).*cos(tau)*AA + AA)./(beta*cos(theta)*(1 + (sin(tau)*gamma - cos(tau)).*exp(tau*gamma)));
z0=-cot(theta) - (beta*(exp(-lambda*tau/beta) - 1) - AA*exp(tau*gamma).*sin(tau))./(beta*sin(theta)*cos(theta).*phi_2(tau));


figure
plot(tau/pi,y0,'r-','displayname','$y_0$')
hold on
plot(tau/pi,z0,'k-','displayname','$z_0$')
xlabel('$\tau(\pi)$','interpreter','latex')
legend('interpreter','latex')

% ratio of the poincare
r=(AA*exp(tau*(2*beta*gamma + lambda)/beta)+...
    ((beta*gamma^2 - (AA + lambda)*gamma + lambda*gamma + beta)*sin(tau) - cos(tau)*AA).*exp((beta*gamma + lambda)*tau/beta)+ ...
    ((-beta*gamma^2 - gamma*lambda - beta)*sin(tau) + lambda*cos(tau)).*exp(tau*gamma) - lambda)./(phi_2(tau)*(AA + lambda));

figure(11)
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
figure(11)
plot(tau/pi,rz,'p-','displayname','$ratio_Z(\tau)$')
legend('interpreter','latex')
%
figure
plot(tau/pi,r,'s-','displayname','$ratio_r(\tau)$')
legend('interpreter','latex')
%
figure
plot(tau/pi,rz,'p-','displayname','$ratio_Z(\tau)$')
legend('interpreter','latex')

% detect the pole of the denom of r_Z
nr_Z = N1+N2+N3;
dr_Z = D1+D2+D3;
figure
plot(tau/pi,dr_Z,'p-','displayname','$ratio_Z(\tau)$')
hold on
plot(tau/pi,nr_Z,'s-','displayname','$ratio_Z(\tau)$')
legend('interpreter','latex')

figure
plot(tau/pi,Asyp1,'p-','displayname','$Asymp1(\tau)$')
hold on
plot(tau/pi,Asyp2,'p-','displayname','$Asymp2(\tau)$')
legend('interpreter','latex')


