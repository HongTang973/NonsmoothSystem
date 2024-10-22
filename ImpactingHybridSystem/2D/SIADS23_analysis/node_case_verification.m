clc
close all
clear
%
lambda_1 = -0.1;
lambda_2 = -0.2;
r        = 2;
det_ = @(t) -r.*(lambda_1 - lambda_2).*exp(t.*(lambda_1 + lambda_2)) + (r.*lambda_1 + lambda_2).*exp(t.*lambda_1) + (-r.*lambda_2 - lambda_1).*exp(t.*lambda_2) + lambda_1 - lambda_2;
%
t_range  = 0:0.01:100;
figure
plot(t_range, det_(t_range))

A        = [lambda_1 + lambda_2,1; -lambda_1*lambda_2, 0];
reverse_R = @(t) (exp(t.*(lambda_1 + lambda_2)).*lambda_1 - exp(t.*(lambda_1 + lambda_2)).*lambda_2 - exp(t.*lambda_1).*lambda_1 + lambda_2.*exp(t.*lambda_2))./(exp(t.*lambda_1).*lambda_2 - lambda_1.*exp(t.*lambda_2) + lambda_1 - lambda_2);
mp = @(t) exp(t.*(lambda_1 + lambda_2))./reverse_R(t).^2;
V  = @(t) (-exp(t.*lambda_1).*lambda_2 + lambda_1.*exp(t.*lambda_2) - lambda_1 + lambda_2)./(exp(t.*lambda_2) - exp(t.*lambda_1));
equi_type = 1;
T        = 10;
% B        =  [0; (1+0.2*1/reverse_R(T))];
B        =  [0; 2];
C        =  [1 0];
R        = eye(2) - B*C*A;
int_2    =  equi_type*[1; 0.98*V(T)-(lambda_1 + lambda_2)];
t_2      = 100;
fs       = 20;

[tout1,yout1,~,teout1,~,~]=...
    Single_DS_Impacting_Hybrid_system_integration(A,R,C,int_2,[0 t_2],fs,equi_type);
figure
plot(tout1,yout1(:,1))

%% 
[v,d] = eig(A);
%
scale  = 2;
figure
plot(v(1,1)*[-1 1]*scale, v(2,1)*[-1 1]*scale, 'k-')
hold on
plot(v(1,2)*[-1 1]*scale, v(2,2)*[-1 1]*scale, 'k-')
plot(yout1(:,1), yout1(:,2))