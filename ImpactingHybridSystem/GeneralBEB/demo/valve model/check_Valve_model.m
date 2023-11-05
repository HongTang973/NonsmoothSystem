% show existence of the LCO
% define the system
clc
clear
close all
% define the linear part

% A_1 = [0,1,0,0,0;
%     37.8653609321384,-0.0499999999989953,2.72979999999201,2.72979999999201,0;
%     -0.00897318244605328,0,-0.000399373906615191,-0.000399373906615191,-1.96243910136930e-05;
%     0.0126899963128540,0,0.000564799995196163,0.000564799995196163,0.0371673910691084;
%     0.108810124583376,-646.644082183967,-0.0143649622907152,-706.385169859214,-1.06913169590479];
A_1 =[0    1.0000         0         0         0;
   37.7917   -0.0500    2.7298    2.7298         0;
   -0.0249         0   -0.0011   -0.0011   -0.0001;
    0.0352         0    0.0016    0.0016    0.4747;
   -0.1290 -280.6388   -0.0058 -221.1266   -5.9395];
[V1,D1]=eig(A_1)
% define initial condition: grid on the positive set
%
state0 = [0.489582485371106;0;5.50000000000000;0;0];
F_ls= @(t,y) A_1*(y-state0);
%
% define the reset map

pr=0.2; %coef of restitution
pz = 10;
C  = [1 0 0 0 0];
vx_= @(t,y)  C*A_1*(y-state0); 
B= -(1+pr)*[0;  1; 0 ; 0; 0];
MW = B*C*A_1;
sign_V = @(Y) sign(C*A_1*Y);
% 
P = eye(size(A_1,1)) + MW;
EA = @(T) P*V1*diag(exp(diag(D1)*T))*inv(V1);
%
A_s = (eye(5)-(B*C*A_1)/(C*A_1*B))*A_1;
eig(A_s)
% for different Evolution time 
 T = 0.0:0.002:2*pi;
%  T = 0.0:0.01:0.1;
 MAX = zeros(1,length(T));
 F_1 = zeros(1,length(T));
 V_sign = zeros(1,length(T));
 LOCI= zeros(5,length(T));
 equi_type =1;
 for i=1:length(T)
     [V_sign(i),LOCI(:,i),MAX(i),F_1(i)] = solve_classification(T(i),EA,sign_V);
%      [sign_V,LOCI,Max,vector,F_1(i),COND] = LCO_Det_search(T(i),P,A_1,C,equi_type);
 end
figure
 plot(T/pi,F_1,'r.','linewidth',1.2,'displayname','F1')
 hold on
%  plot(beta*T/pi,V_sign,'-','linewidth',1.2,'displayname','Sign of Velocity')
 plot([0 2],[0 0],'b-','displayname','zero line')
%  plot(beta*T/pi,MAX,'k-','linewidth',1.0,'displayname','MAX-1')
%  plot(0.986964041540607,0,'g*')
% beta =1;
%  plot(beta*T/pi,LOCI(1,:),'y--','linewidth',1.4,'displayname','EIG1-1')
%  plot(beta*T/pi,LOCI(2,:),'g--','linewidth',1.4,'displayname','EIG1-2')
%  plot(beta*T/pi,LOCI(3,:),'r--','linewidth',1.4,'displayname','EIG1-3')
 legend
%  ylim([-2 2])
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
M = EA(T)-eye(5);
PA = EA(T);
Ans = -(M(2:5,2:5))\M(2:5,1);
F_1 = M(1,:)*[1;Ans];
[JV,JD]=eig(PA);
[MaxEig,ind] = max(abs(diag(JD)));
v=JV(:,ind);
v = v/v(1);
Max = MaxEig-1;
LOCI= sort(abs(diag(JD)))-1;
sign =sign_V(v);
end



