% singular matrix system

clc
clear
close all
par = [0;0.00730000000000000;-1;2.01707011290029;1;0.100000000000000];

lambda_1 = par(1);

lambda_2 = par(2);

lambda_3 = par(3);

a1 = lambda_1 + lambda_2 + lambda_3;
a2 = -(lambda_1 * lambda_2 + lambda_2*lambda_3 + lambda_1*lambda_3);
a3 = lambda_1 * lambda_2 * lambda_3;
%
b2 = par(4);
b3 =par(5);

A = [a1 1 0;a2 0 1; a3 0 0];


B = [0; b2; b3];


C = [1,0,0];

R = eye(length(C)) - B*C*A;

T = par(6);

equi_type = -1;
fs = 10;

[~,~,~,LCO,~] = LCO_Det_search(T,R,A,C,-equi_type);
InitCond = LCO;
InitCond = [1;1;0.0073];
tol.tol_v = 1e-3*abs(C*A*LCO);
[tout,yout,yeout0,teout,yeout,ieout]=...
SDS_IHS_Int(A,R,C,sign(C*A*LCO)*InitCond,[0 1000],fs,sign(C*A*LCO),tol);

figure
plot(yout(:,1),yout(:,2),'k-','linewidth',1.2)

figure
plot(tout,yout(:,1),'k-','linewidth',1.2)

figure
plot(tout,yout(:,2),'k-','linewidth',1.2)