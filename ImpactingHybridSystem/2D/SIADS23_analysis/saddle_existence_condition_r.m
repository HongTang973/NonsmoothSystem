clc
clear
close all

lambda_1 = 0.1;
lambda_2 = -0.2;

reverse_R = @(t) (exp(t.*(lambda_1 + lambda_2)).*lambda_1 - exp(t.*(lambda_1 + lambda_2)).*lambda_2 - exp(t.*lambda_1).*lambda_1 + lambda_2.*exp(t.*lambda_2))./(exp(t.*lambda_1).*lambda_2 - lambda_1.*exp(t.*lambda_2) + lambda_1 - lambda_2);
mp        = @(t) exp(t.*(lambda_1 + lambda_2))./reverse_R(t).^2;
V         = @(t) (-exp(t.*lambda_1).*lambda_2 + lambda_1.*exp(t.*lambda_2) - lambda_1 + lambda_2)./(exp(t.*lambda_2) - exp(t.*lambda_1));
t_range   = 0:0.01:100;

figure
plot(t_range, reverse_R(t_range),'k')
hold on
% plot(t_range, mp(t_range),'r--')
% plot(t_range, V(t_range),'b--')
