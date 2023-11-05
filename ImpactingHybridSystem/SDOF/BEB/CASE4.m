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

public.plot.addmissible_ob.color = [0.2 ,0.2 , 0.2]; %black
public.plot.addmissible_ob.linewidth = 1;
public.plot.addmissible_ob.style= '-';
%
public.plot.S_po.color = [0 ,0 , 1]; % blue
public.plot.S_po.linewidth = 1.8;
public.plot.S_po.style = '-';
%
public.plot.US_po.linewidth = 2;
public.plot.US_po.color = [1,0,0];
public.plot.US_po.style = '--';

%
public.plot.gca.linewidth = 3;
public.plot.gca.fontname = 'Times New Roman';
public.plot.gca.fontsize = 14;
%
HX = @(C,Y,mu) C*Y - mu;
HV = @(C,A,Y) C*A*Y;
public.function.HX = HX;
public.function.HV = HV;
%
BEB2D.public =public;

%% -----------CASE 4: From admissible focus to pseudo focus --------------%
% stable sdmissible equilibrium when mu =-1;
% stable LCO when mu =1;

%------ define the system ------%
CASE = [];
CASE.num = 1;
CASE.description = 'AF2PF_po_S';

a =0.5 ;b=-1;
A = [0,1;b,a];
r= 0.4;
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



% DNS for the  ------------ mu=-1 ----------
equi_type =   -1; 

MINUS = [];
MINUS.description = 'addmissible unstable equilibrium and stable LCO';
MINUS.mu = -1;
fig_num =7;

% presicted by technique developed in paper
MINUS.po.starter = equi_type*[1;-7.0237];
MINUS.po.period = 3.44;
MINUS.po.Mono_p = [1;-2.2327];
MINUS.po.Salt_p = [1;0.8931];
MINUS.po.description = 'Stable LCO';
MINUS.mu =  -1;
MINUS.fs =fs;


% post processing of the data
yt=0;
yb=0;
xr=0;
% try some plot in this case
int_2 = MINUS.po.starter;
 % outside LCO
int_2(2) = int_2(2) + 2;
t_2 = 30*pi;
[tout1,yout1,~,teout1,~,~]=...
    Single_DS_Impacting_Hybrid_system_integration(A,R,C,int_2,[0 t_2],fs,equi_type);
%
MINUS.init2 = int_2;
MINUS.tspan2 = t_2;


[ytop,ybot,xright,CASE,BEB2D]=trunc_plot(tout1,yout1,teout1,floor(t_2),CASE,MINUS,fig_num,BEB2D);

yt = max([yt,ytop]);
yb = min([yb,ybot]);
xr = max([xr,xright]);
%
int_1 = MINUS.po.starter;
int_1(2) = int_1(2) -2;
t_1 = 30*pi;
[tout1,yout1,~,teout1,~,~]=...
    Single_DS_Impacting_Hybrid_system_integration(A,R,C,int_1,[0 t_1],fs,equi_type);
%
MINUS.init1 = int_1;
MINUS.tspan1 = t_1;

[ytop,ybot,xright,CASE,BEB2D]=trunc_plot(tout1,yout1,teout1,floor(t_1),CASE,MINUS,fig_num,BEB2D);

%
int_1 = MINUS.po.starter;
int_1(2) = int_1(2) -2;
t_1 = 30*pi;
[tout1,yout1,~,teout1,~,~]=...
    Single_DS_Impacting_Hybrid_system_integration(A,R,C,int_1,[0 t_1],fs,equi_type);
%
MINUS.init1 = int_1;
MINUS.tspan1 = t_1;

[ytop,ybot,xright,CASE,BEB2D]=trunc_plot(tout1,yout1,teout1,floor(t_1),CASE,MINUS,fig_num,BEB2D);

yt = max([yt,ytop]);
yb = min([yb,ybot]);
xr = max([xr,xright]);

% EXACT　LCO
int_0 = MINUS.po.starter;
[~,PO_S,~,~,~,~]=...
    Single_DS_Impacting_Hybrid_system_integration(A,R,C,int_0,[0 MINUS.po.period],fs,equi_type);
MINUS.po.seg = PO_S;

[X,Y] = Observing_filter(public,CASE,PO_S',MINUS.mu);
Y(X<0)=[];
X(X<0)=[];
 % ---- Plot---
 
 MINUS.xlim = [-1,18];
 MINUS.ylim = [-25 , 15];
 
figure(fig_num)
plot(feval(public.function.HX,C,[0;0],MINUS.mu),0,'ro','LineWidth', 1.5,'markerfacecolor',[1,1,1])
hold on

plot([0,0],25*[-1 1],'--','color',[0,0,0])
plot(X,Y,'linestyle',public.plot.S_po.style,'linewidth',public.plot.S_po.linewidth,'color',public.plot.S_po.color)

% xlabel('$H(x)$','interpreter','latex')
% ylabel('$\frac{\rm{d}}{\rm{d} t}H(x)$','interpreter','latex')

xlim(MINUS.xlim)
ylim(MINUS.ylim)
set(gca,'fontsize',public.plot.gca.fontsize,...
        'fontname',public.plot.gca.fontname,...
        'linewidth',public.plot.gca.linewidth)      

% ------------------ mu =1 ----------------------------%
equi_type =   1; 
POS = [];

POS.description = 'stable pseudo equilibrium';
POS.mu =  1;
POS.fs =fs;
fig_num = 8;
% post processing of the data

yt=0;
yb=0;
xr=0;
% try some plot in this case
int_1 = zeros(2,1);
int_1(1) = equi_type;
int_1(2) = 4;
t_1 = 15*pi;
[tout1,yout1,~,teout1,~,~]=...
    Single_DS_Impacting_Hybrid_system_integration(A,R,C,int_1,[0 t_1],fs,equi_type);
figure(1)
plot(yout1(:,1),yout1(:,2),'k--','linewidth',0.8)
%
POS.init = int_1;
POS.tspan = t_1;
POS.fs =fs;

[ytop,ybot,xright,CASE,BEB2D]=trunc_plot(tout1,yout1,teout1,floor(t_1),CASE,POS,fig_num,BEB2D);
yt = max([yt,ytop]);
yb = min([yb,ybot]);
xr = max([xr,xright])-equi_type+0.5;

figure(fig_num)
plot(feval(public.function.HX,C,[1;0],POS.mu),0,'ro','LineWidth', 1.5,'markerfacecolor',[0,0,0])
hold on
plot([0,0],25*[-1 1],'--','color',[0,0,0])
% xlabel('$H(x)$','interpreter','latex')
% ylabel('$\frac{\rm{d}}{\rm{d} t}H(x)$','interpreter','latex')
xlim([-1.1,6])
ylim([-10 10])
set(gca,'fontsize',public.plot.gca.fontsize,...
        'fontname',public.plot.gca.fontname,...
        'linewidth',public.plot.gca.linewidth)
