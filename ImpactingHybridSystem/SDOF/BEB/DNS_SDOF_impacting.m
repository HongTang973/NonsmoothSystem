% Detecting the existence of LCO when boundary equilibrium transferring to
% virtual Equilibrium for general n dimension impacting hybrid system 
% limitations: single impact; single discontinuity surface

% KEY INPUTs:
% A : the linear region system's matrix 
% R : the reset map matrix, R(x) = R*x;
% C : the selecting vector to define the discontinuity surface C*x=0

clc
clear
close all
% define the linear part

% xi = 1.05; r = 1.5;
% A  =  [0 1;-1 -2*xi]; 
BEB2D =[];

public.plot.addmissible_ob.color = [0 ,0 , 0]; %black
public.plot.addmissible_ob.linewidth = 1;
public.plot.addmissible_ob.style= 'solid';
%
public.plot.S_po.color = [0 ,0 , 1]; % blue
public.plot.S_po.linewidth = 1.8;
public.plot.S_po.style = 'solid';
%
public.plot.gca.linewidth = 1.5;
public.plot.gca.fontname = 'Times New Roman';
public.plot.gca.fontsize = 12;
%
HX = @(C,Y,mu) C*Y - mu;
HV = @(C,A,Y) C*A*Y;
public.function.HX = HX;
public.function.HV = HV;
%
BEB2D.public =public;

%% -----------CASE 1: From admissible focus to pseudo focus --------------%
% stable sdmissible equilibrium when mu =-1;
% stable LCO when mu =1;

%------ define the system ------%
CASE = [];
CASE.num = 1;
CASE.description = 'AF2PF_po_S';

a =-1 ;b=-1;
A = [0,1;b,a];
r= 1.5;
MW = [0 0;0 -(1+r)];
R  = eye(size(A,1)) + MW;
C  = [1,0];
%
CASE.C=C;
CASE.A=A;
CASE.P=R;
% define the useful function
[V1,D1]=eig(A);
CASE.eigSOL.V1 = V1;
CASE.eigSOL.D1 = D1;

EA = @(T) real(V1*diag(exp(diag(D1)*T))*inv(V1));
sign_V = @(Y) sign(C*A*Y);

% determine the sampling frequency 
Omega = max(abs(diag(D1)));
fs = 100*ceil(2*(Omega/2/pi));
t0 =150;
tspan = [0 t0];

% DNS   
% if equi_type == 'admissible'
%         pre = -1;
%  elseif equi_type == 'pseudo'
%         pre = 1;
%  end

% classification condition

CASE.tr = sign(sum(diag(D1)));
CASE.rexp = r*exp(real(D1(1,1))*pi/imag(D1(1,1)));
fprintf('The condition trace(A) sign:\n %d \n',CASE.tr);
fprintf('The condition r*exp(a*pi/w):\n %f \n',CASE.rexp);



% DNS for the 
equi_type =   -1; 

MINUS = [];
MINUS.description = 'addmissible stable equilibrium';
MINUS.mu = -1;
fig_num =1;

% post processing of the data

yt=0;
yb=0;
xr=0;
% try some plot in this case
int_1 = zeros(2,1);
int_1(1) = equi_type;
int_1(2) = 15;
t_1 = 10*pi;
[tout1,yout1,~,teout1,~,~]=...
    Single_DS_Impacting_Hybrid_system_integration(A,R,C,int_1,[0 t_1],fs,equi_type);
% figure(1)
% plot(yout1(:,1),yout1(:,2),'k--','linewidth',0.8)
%
MINUS.init = int_1;
MINUS.tspan = t_1;
MINUS.fs =fs;

[ytop,ybot,xright,CASE,BEB2D]=trunc_plot(tout1,yout1,teout1,floor(t_1),CASE,MINUS,fig_num,BEB2D);
yt = max([yt,ytop]);
yb = min([yb,ybot]);
xr = max([xr,xright])-equi_type+0.5;

figure(fig_num)
plot(feval(public.function.HX,C,[0;0],MINUS.mu),0,'ro','LineWidth', 1.5,'markerfacecolor',[0,0,0])
hold on
plot([0,0],10*[-1 1],'--','color',[0,0,0])
xlabel('$H(x)$','interpreter','latex')
ylabel('$\frac{\rm{d}}{\rm{d} t}H(x)$','interpreter','latex')
xlim([-1.1,xr])
set(gca,'fontsize',public.plot.gca.fontsize,...
        'fontname',public.plot.gca.fontname,...
        'linewidth',public.plot.gca.linewidth)

% ------------------ mu =-1 ----------------------------%
equi_type =   1; 
POS = [];

POS.description = 'stable periodic orbit';
% presicted by technique developed in paper
POS.po.starter = [1;0.7950];
POS.po.period = 0.3734 * pi;
POS.po.Mono_p = [1;-0.4641];
POS.po.Salt_p = [1;0.6962];
POS.po.description = 'Stable LCO';
POS.mu =  1;
POS.fs =fs;

% post processing of the data
fig_num =2;
yt=0;
yb=0;
xr=0;
% try some plot in this case
int_2 = POS.po.starter;
 % outside LCO
int_2(2) = int_2(2) + 0.5;
t_2 = 4*pi;
[tout1,yout1,~,teout1,~,~]=...
    Single_DS_Impacting_Hybrid_system_integration(A,R,C,int_2,[0 t_2],fs,equi_type);
%
POS.init2 = int_2;
POS.tspan2 = t_2;


[ytop,ybot,xright,CASE,BEB2D]=trunc_plot(tout1,yout1,teout1,floor(t_2),CASE,POS,fig_num,BEB2D);

yt = max([yt,ytop]);
yb = min([yb,ybot]);
xr = max([xr,xright]);
% inside LCO
int_1 = POS.po.starter;
int_1(2) = int_1(2) -0.5;
t_1 = 5*pi;
[tout1,yout1,~,teout1,~,~]=...
    Single_DS_Impacting_Hybrid_system_integration(A,R,C,int_1,[0 t_1],fs,equi_type);
%
POS.init1 = int_1;
POS.tspan1 = t_1;


[ytop,ybot,xright,CASE,BEB2D]=trunc_plot(tout1,yout1,teout1,floor(t_1),CASE,POS,fig_num,BEB2D);

yt = max([yt,ytop]);
yb = min([yb,ybot]);
xr = max([xr,xright]);
% EXACT　LCO
int_0 = POS.po.starter;
[~,PO_S,~,~,~,~]=...
    Single_DS_Impacting_Hybrid_system_integration(A,R,C,int_0,[0 POS.po.period],fs,equi_type);
POS.po.seg = PO_S;

[X,Y] = Observing_filter(public,CASE,PO_S',POS.mu);
 % ---- Plot---
 
 POS.xlim = [-0.2,0.4];
 POS.ylim = [-2 , 2];
 
figure(fig_num)
plot(feval(public.function.HX,C,[MINUS.mu;0],MINUS.mu),0,'ro','LineWidth', 1.5,'markerfacecolor',[1,1,1])
hold on
plot([0,0],10*[-1 1],'--','color',[0,0,0])
plot(X,Y,'-','linewidth',public.plot.S_po.linewidth,'color',public.plot.S_po.color)

xlabel('$H(x)$','interpreter','latex')
ylabel('$\frac{\rm{d}}{\rm{d} t}H(x)$','interpreter','latex')
xlim(POS.xlim)
ylim(POS.ylim)
set(gca,'fontsize',public.plot.gca.fontsize,...
        'fontname',public.plot.gca.fontname,...
        'linewidth',public.plot.gca.linewidth)
    
%% CASE2:


    
    
    
    
%  plot(equi_type,0,'ro','markerfacecolor',[0,0,0])




% plot for the node/sadle asmptotic line 
scalar = -5:0.01:8;
line_1 = V1(:,1)*scalar;
line_2 = V1(:,2)*scalar;
right_1 = find(line_1(1,:)>equi_type);
right_2 = find(line_2(1,:)>equi_type);
line_1 = line_1(:,right_1);
line_2 = line_2(:,right_2);

% figure
% plot(line_1(1,:),line_1(2,:),'-','color',[0,0,0],'linewidth',1.4)
% hold on
% plot(line_2(1,:),line_2(2,:),'-','color',[0,0,0],'linewidth',1.4)

%
% InitCond = [equi_type ]';
% [tout,yout,yeout0,teout,yeout,ieout]=...
%     Single_DS_Impacting_Hybrid_system_integration(A,R,C,InitCond,tspan,fs,equi_type);
% Period = teout(2:end)-teout(1:end-1);
% %
% figure
% plot(tout,yout(:,C>0),'r-','linewidth',1.4)
% title(['LCO ',num2str(i),'''s stability'])
% set(gca,'fontname','times new roman','fontsize',12)
% xlabel('t/s')
% figure
% plot(yout(:,1),yout(:,2),'k-','linewidth',1.2)
% title(['LCO ',num2str(i),'''s phase portrait'])
% set(gca,'fontname','times new roman','fontsize',12)
% grid on


% [ytop,ybot,xright]=trunc_plot(tout,yout,teout,'blue',20,fs,fig_num,1.8);
% yt = max([yt,ytop]);
% yb = min([yb,ybot]);
% xr = max([xr,xright]);


% auxilary line to show the basins
int_1 = zeros(2,1);
int_1(1) = equi_type;
int_1(2) = 5.9;
t_1 = 100*pi;
[tout1,yout1,~,teout1,~,~]=...
    Single_DS_Impacting_Hybrid_system_integration(A,R,C,int_1,[0 t_1],fs,equi_type);
[ytop,ybot,xright]=trunc_plot(tout1,yout1,teout1,[0.5,0.5,0.5],t_1,fs,fig_num,0.8);
yt = max([yt,ytop]);
yb = min([yb,ybot]);
xr = max([xr,xright]);
%
int_2 = zeros(2,1);
int_2(1) = 0*equi_type;
int_2(2) = -1.;
t_2 = 0.3*pi;
[tout2,yout2,~,teout2,~,~]=...
    Single_DS_Impacting_Hybrid_system_integration(A,R,C,int_2,[0 t_2],fs,equi_type);
 [ytop,ybot,xright]=trunc_plot(tout2,yout2,teout2,[0.5,0.5,0.5],t_2,fs,fig_num,0.8);
yt = max([yt,ytop]);
yb = min([yb,ybot]);
xr = max([xr,xright]);
keyboard
% plot(equi_type-0.075*abs(xright-equi_type),0,'ro')
% plot([equi_type,equi_type],500*[-1 1],'--','color',[0,0,0])
% xlim([equi_type-0.15*abs(xright-equi_type),xright+0.1*abs(xright-equi_type)])
% ylim([ybot-0.15*abs(ytop-ybot),ytop+0.15*(ytop-ybot)])
%
int_3 = [-1.1*equi_type, -1]';
% t_3 = 0.9503*pi; tt = 0:1/fs:t_3;
t_3 = 0.3*pi;
[tout3,yout3,~,teout3,~,~]=...
    Single_DS_Impacting_Hybrid_system_integration(A,R,C,int_3,[0 t_3],fs,equi_type);
% % ayout=[];
% % for kk = 1:length(tt)    
% % ayout = [ayout,EA(tt(kk))*int_3];
% % end
figure(1)
plot(yout3(:,1),yout3(:,2),'k--','linewidth',0.8)
 hold on
 plot(ayout(1,:),ayout(2,:),'-b','linewidth',1.8)
figure(fig_num)
% [ytop,ybot,xright]=trunc_plot(tout,yout,teout,[0.5,0.5,0.5],t0,fs,fig_num,0.8);
%  plot(equi_type,0,'ro','markerfacecolor',[0,0,0])
plot(0,0,'ro')
hold on
plot([equi_type,equi_type],10*[-1 1],'--','color',[0,0,0])
set(gca,'fontname','times new roman','xticklabel',[],'yticklabel',[])




