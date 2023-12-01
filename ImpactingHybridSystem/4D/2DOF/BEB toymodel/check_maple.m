% check the fidelity of the results
clc
clear
close all
%
alpha = -0.1;
beta = 0.2;
y0 = 0.53;
t=0:0.01:100;
x =  exp(alpha.*t).*cos(beta.*t) + exp(alpha.*t).*sin(beta.*t).*y0-1;
y = -exp(alpha.*t).*sin(beta.*t) + exp(alpha.*t).*cos(beta.*t).*y0;

figure
plot(t,x)

figure
plot(t,y)