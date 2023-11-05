% show existence of the LCO
% define the system
clc
clear
close all
% define the linear part
alpha = -0.1;
beta  = 0.2;
lambda = -0.5;
gamma = (alpha-lambda)/beta;
A_1=[alpha beta 0; -beta alpha 0;0 0 lambda];
theta = pi/6;
P = [cos(theta),0,-sin(theta);0,1,0;sin(theta),0,cos(theta)];
A_1 = inv(P)*A_1*P;
[V1,D1]=eig(A_1);
% define initial condition: grid on the positive set
%
state0 = [-1 0 0]';
F_ls= @(t,y) A_1*(y-state0);
%
% define the reset map
a11=A_1(1,1);a12=A_1(1,2);a13=A_1(1,3);
N = sqrt(a12^2+a13^2);
V_l3 = [0;a13;-a12]/N;
P_l3 = [0;a12; a13]/N;
pr=0.8; %coef of restitution
pz = 10;
vx_= @(t,y)  [1 0 0]*A_1*(y-state0); 
Wx_= (-(1+pr)*[0 a12 a13]'+pz*[0 a13 -a12]')./(a12^2+a13^2);
MW =  -(1+pr)*P_l3*[a11 a12 a13]/N+pz*V_l3*[a11 a12 a13]/N;
sign_V = @(Y) sign([1 0 0]*A_1*Y);
% 
P = eye(size(A_1,1)) + MW;
EA = @(T) P*V1*diag(exp(diag(D1)*T))*inv(V1) ;
%
% for different Evolution time 
 T = 0.0*pi/beta:0.002:2*pi/beta;
%  T = 0.0:0.01:0.1;
 MAX = zeros(1,length(T));
 F_1 = zeros(1,length(T));
 V_sign = zeros(1,length(T));
 LOCI= zeros(3,length(T));
 for i=1:length(T)
     [V_sign(i),LOCI(:,i),MAX(i),F_1(i)] = solve_classification(T(i),EA,sign_V);
 end
figure
 plot(beta*T/pi,F_1,'r-','linewidth',1.2,'displayname','F1')
 hold on
%  plot(beta*T/pi,V_sign,'-','linewidth',1.2,'displayname','Sign of Velocity')
 plot([0 2],[0 0],'b-','displayname','zero line')
%  plot(beta*T/pi,MAX,'k-','linewidth',1.0,'displayname','MAX-1')
%  plot(0.986964041540607,0,'g*')
 plot(beta*T/pi,LOCI(1,:),'y--','linewidth',1.4,'displayname','EIG1-1')
 plot(beta*T/pi,LOCI(2,:),'g--','linewidth',1.4,'displayname','EIG1-2')
 plot(beta*T/pi,LOCI(3,:),'r--','linewidth',1.4,'displayname','EIG1-3')
 legend
 ylim([-2 2])
 xlabel('\tau /\pi')
 ylabel('value')
 set(gca,'fontname','Times New Roman')
 
 % verification the slope 
%  tau = 0:0.01:2*pi;
%  flow_dir = -(sin(theta)*pz*gamma + pr + 1)/((1 + pr)*sin(theta)*gamma - pz);
%  slope = -sin(theta)*(cos(tau) - exp(-tau*gamma))./sin(tau) -flow_dir;
%  figure
%  plot(tau/pi,slope);
%  plot([0 2],[0 0])
%  hold on
%  plot(beta*T/pi,MAX,'k-','linewidth',1.2,'displayname','MAX-1')
%  ylim([-2 2])
% [~,Max,F_1] = solve_classification(1.29365*pi/beta,EA)

function [sign,LOCI,Max,F_1] = solve_classification(T,EA,sign_V)
M = EA(T)-eye(3);
PA = EA(T);
Ans = -(M(2:3,2:3))\M(2:3,1);
F_1 = M(1,:)*[1;Ans];
[JV,JD]=eig(PA);
[MaxEig,ind] = max(abs(diag(JD)));
v=JV(:,ind);
v = v/v(1);
Max = MaxEig-1;
LOCI= sort(abs(diag(JD)))-1;
sign =sign_V(v);
end



