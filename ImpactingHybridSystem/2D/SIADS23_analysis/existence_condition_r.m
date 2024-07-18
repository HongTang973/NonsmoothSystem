clc
clear
close all
mu= 0.1;
x_0=-1;
y_0=0.1;
x=@(tau)  exp(-mu.*tau).*(x_0*cos(tau)+(y_0*sqrt(1+mu^2)+mu*x_0).*sin(tau));
y=@(tau)  exp(-mu*tau).*(y_0*cos(tau)-(x_0*sqrt(1+mu^2)+mu*y_0).*sin(tau));
F=@(tau) 1-exp(mu*tau).*(cos(tau)-mu*sin(tau));
%> define the auxilary function
phi=    @(tau,mu) 1-exp(mu*tau).*(cos(tau)-mu*sin(tau));
%> 
r  = @(tau, mu)  phi(tau, mu)./phi(tau, -mu).*exp(-2*mu*tau);
mp = @(tau, mu)  phi(tau, mu).^2./phi(tau, -mu).^2.*exp(-2*mu*tau);
sqrt_mp = @(tau,mu)  phi(tau, mu)./phi(tau, -mu).*exp(-mu*tau);
t_range = -2*pi:0.01:2*pi;
figure
plot(t_range, r(t_range,mu), 'k-')
hold on 
plot(t_range, mp(t_range,mu), 'r-')
plot(t_range, sqrt_mp(t_range,mu), 'g-')
plot(pi, r(pi,mu), 'bs')
plot(-pi, r(-pi,mu), 'rs')
plot(0, 1, 'rs')
ylim([-25 25])
figure
plot(t_range,F(t_range))


dd_sqrt_mp = @(T,mu) (-cos(T).^2 + 4.*sin(T).*mu.*(mu.^2 + 1).*cos(T) - mu.^4.*sin(T).^2 - 2.*mu.^2 + 2).*exp(2.*T.*mu) + (mu.*sin(T) - cos(T)).*exp(3.*T.*mu) + (mu.^2.*cos(T).^3 + cos(T).^2.*sin(T).*mu.^3 + (-mu.^4.*sin(T).^2 + 3.*mu.^2 + 1).*cos(T) - sin(T).*mu.*(mu.^4.*sin(T).^2 + 5.*mu.^2 + 7)).*exp(T.*mu) + (-mu.^2 + 1).*cos(T).^2 + 2.*mu.*cos(T).*sin(T) - mu^2 - 2;
d_half_mp  = @(T,mu) (cos(T).*mu + sin(T)).*exp(2.*T.*mu) + (-mu.^3.*sin(T).^2 + mu.*cos(T).^2 - 3.*mu).*exp(T.*mu) + cos(T).*mu - sin(T);
dd_half_mp = @(T,mu) (2.*mu.^2.*cos(T) + mu.*sin(T) + cos(T)).*exp(2.*T.*mu) - mu.*(-mu.*cos(T).^2 + (2.*mu.^2 + 2).*sin(T).*cos(T) + mu.^3.*sin(T).^2 + 3.*mu).*exp(T.*mu) - mu.*sin(T) - cos(T);
figure
hold on 
plot(t_range, sqrt_mp(t_range,mu), 'g-')
plot(t_range, ones(size(t_range)),'b-')
plot(t_range, d_half_mp(t_range,mu)./phi(-t_range,mu).^2./exp(t_range.*mu).^2, 'r-')
% plot(t_range, dd_half_mp(t_range,mu),'r-')
% plot(t_range, -dd_sqrt_mp(t_range,mu))
% plot(t_range, -dd_sqrt_mp(t_range,mu)./phi(-t_range,mu).^3./exp(t_range.*mu).^3)
plot(t_range, 0*t_range,'b-')
plot([0.5*pi 0.5*pi],[-15 15], 'k-')
plot(-[0.5*pi 0.5*pi],[-15 15], 'k-')
ylim([-1 2])
xlim([0 1.5*pi])

%% 
% test = @(T,mu) (-cos(T).*mu - sin(T)).*exp(2.*T.*mu) + ((mu.^2 + 2).*sin(T).^2 + cos(T).^2 + 1).*mu.*exp(T.*mu) - cos(T).*mu + sin(T);
% figure
% plot(t_range, test(t_range,mu))
test = @(T,mu) (sin(T).*mu.^2 - mu + sin(T)).*exp(T.*mu) + cos(T).*mu - sin(T);
plot(t_range, test(t_range,mu),'b-')