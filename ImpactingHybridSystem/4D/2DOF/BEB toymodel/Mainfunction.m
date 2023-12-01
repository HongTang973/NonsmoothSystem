% x_0: incoming point
% x_00: outcoming point
% x_p : point on the poincare section
% T: the period of this orbit
% F_ls: function handle of the flow
% RMap: the reset map when the flow hit the boundary
% EventFun: the function to detect the event of hitting the boundary
% PhaseInd: 1x3 vector, with the first element  specifying the obsevor 
% last two variables to plot out phase 
function Mainfunction
% define the linear part
alpha = -0.1;
beta  = 0.2;
lambda = -0.5;
A_1=[alpha beta 0; -beta alpha 0;0 0 lambda];
theta = pi/6;
P = [cos(theta),0,-sin(theta);0,1,0;sin(theta),0,cos(theta)];
A_1 = inv(P)*A_1*P;
% define initial condition: grid on the positive set
%
state0 = [-1 0 0]';
%
F_ls= @(t,y) A_1*(y-state0);
%
% define the reset map
a11=A_1(1,1);a12=A_1(1,2);a13=A_1(1,3);
N = sqrt(a12^2+a13^2);
V_l3 = [0;a13;-a12]/N;
P_l3 = [0;a12; a13]/N;
pr=0.8; %coef of restitution
pz = 8;

RMap = @(y)  y +(-(1+pr)*P_l3*A_1(1,:)/N+pz*V_l3*A_1(1,:)/N)*(y-state0);
% 
PlotInd = [1 2 3 1 2];
% %
% T=800;
% tspan=[0 T];
% fs=100;
% mu =1;
% x0=0;
% C=0.1;
% z0=2;
% y0 = ((-alpha + lambda)*cos(theta)^2 + z0*(alpha - lambda)*sin(theta)*cos(theta) + C - lambda)/(beta*cos(theta));
% InitCond=[x0 y0 z0]';
% [~,yout,teout,yeout,yeout0,ieout,amplitude]=PWSC_int(A_1,InitCond,state0,tspan,fs,mu,0);
% Period=beta*(teout(2:end)-teout(1:end-1))/pi;
% %
% T = Period(end)*pi/beta;
% % chose x_p
% cut = yout(end-fs*T:end,:);
% [~,Ind] = max (cut(:,1));
% x_p = cut(Ind,:)';
% load('C:\Users\ib20968\OneDrive - University of Bristol\Codes stall\2DOF\BEB toymodel\CompsedMap_data_set.mat');
load('F:\onedrive\OneDrive - University of Bristol\Codes stall\2DOF\BEB toymodel\CompsedMap_data_set.mat');
EventFun(T,x_p)
fs =100;
clc
close all

% [MainV,LPE]=ComposedMap_FloqueMPs(x_p,T,fs,F_ls,RMap,@EventFun,PlotInd);
% after numerically get the Floque mulitipliers around the LCO, try to find
% them in an analytical way
% construct the Jacobi matrix of the composed map
W =  -(1+pr)*P_l3*[a11 a12 a13]/N+pz*V_l3*[a11 a12 a13]/N;
% T = Period(end)*pi/beta;
[~,~,T_2,x_0,~,~,~]=PWSC_int(A_1,x_p,state0,[0 T],fs,1,1);
% T_2=0;
T_1 = T - T_2;
f_ls= @(t,y) A_1*(y-state0);
[V1,D1]=eig(A_1);
DS=diag(exp(diag(D1)*T_1));
DS0=diag(exp(diag(D1)*T_2));
S  = diag(DS);
S0 = diag(DS0);
KK0=V1*DS0*inv(V1);
KK1=V1*DS*inv(V1);
x_0     = x_0';
x_00    = x_0+W*(x_0-state0);
x_00(1) = 0;
x_0(1)  = 0;
C       = [1 0 0];
P =eye(3)+W;
correction = ((f_ls(0,x_00)-P*f_ls(0,x_0))*C)/(C*f_ls(0,x_0));
Q=P + 1*correction;
J=KK0*Q*KK1
% KK1*(x_00-state0)+state0
% exp(A_1*T)
% inv(W)
[JV,JD] = eig(J)
JV(:,1)'*(x_00-state0)/norm(JV(:,1))/norm(x_00-state0)
% MM = PP*exp(A_1*T);
% [VMM,DMM] = eig(MM)
keyboard

function [value,isterminal,direction] = EventFun(t,y)
% Locate the time when flap hit through boundary in both directions
% and stop integration when isterminal ==1
% Events: switch boundary; zero heave velocity; zero pitch velocity; zero
% flap velocity
value = (y(1));           %  detect switching point
isterminal = 1;                           %  stop the integration
direction  = -1; 
