% check the roots
clc
close all
gamma = -0.1;
phi_1 = @(tau) 1-exp(gamma*tau).*cos(tau);
phi_2 = @(tau) phi_1(tau) -exp(gamma*tau).*sin(tau)/gamma;

tau=0:0.01:2*pi;

p1 = phi_1(tau);
p2 = phi_2(tau);
p3 = 0*tau;

figure
% plot(tau/pi,p1,'r-','displayname','$\varphi_1(\tau)$')
hold on
plot(tau/pi,p3,'b-','displayname','$\rm{zero~line}$')
plot(tau/pi,p2,'k-','displayname','$\varphi_2(\tau)$')
xlabel('$\tau(\pi)$','interpreter','latex')
legend('interpreter','latex')