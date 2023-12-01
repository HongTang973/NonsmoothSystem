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

% define the struct of BEB2D
BEB2D =[];

% define the public contaning the common features
public.plot.addmissible_ob.color = [0 ,0 , 0]; %black
public.plot.addmissible_ob.linewidth = 1;
public.plot.addmissible_ob.style= '-';
%
public.plot.S_po.color = [0 ,0 , 1]; % blue
public.plot.S_po.linewidth = 2;
public.plot.S_po.style = '-';

%
public.plot.US_po.linewidth = 2;
public.plot.US_po.color = [1,0,0];
public.plot.US_po.style = '--';
%

public.plot.gca.linewidth = 3;
public.plot.gca.fontname = 'Times New Roman';
public.plot.gca.fontsize = 16;
public.plot.gca.ws_x = 0.48;
public.plot.gca.ws_y = 0.3354;
public.plot.gca.width = 4.5208;
public.plot.gca.height = 3.5656;




% observing function
HX = @(C,Y,mu) C*Y - mu;
HV = @(C,A,Y) C*A*Y;
public.function.HX = HX;
public.function.HV = HV;

%
BEB2D.public =public;

filepath = [pwd,'\figures'];
if ~exist(filepath, 'dir')
       mkdir(filepath)
end
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

FIG1=figure(fig_num);
% axh1 =axes(FIG1,1);
plot(feval(public.function.HX,C,[0;0],MINUS.mu),0,'ro','LineWidth', 1.5,'markerfacecolor',[0,0,0]);
hold on
plot([0,0],10*[-1 1],'--','color',[0,0,0])
% xlabel('$H(x)$','interpreter','latex')
% ylabel('$\frac{\rm{d}}{\rm{d} t}H(x)$','interpreter','latex')
xlim([-1.1,xr])
set(gca,'fontsize',public.plot.gca.fontsize,...
        'fontname',public.plot.gca.fontname,...
        'linewidth',public.plot.gca.linewidth,...
        'Units','inches',...
        'position',[public.plot.gca.ws_x  public.plot.gca.ws_y ...
     public.plot.gca.width public.plot.gca.height])
    set(FIG1, 'PaperUnits', 'inches','PaperPosition', ...
    [0  0 ...
     public.plot.gca.ws_x+public.plot.gca.width  public.plot.gca.ws_y+public.plot.gca.height],...
    'PaperSize', [public.plot.gca.width public.plot.gca.height]);
    %save 
    name ='2D_BEB_case1_1equi';
    picname1 = strcat(filepath,['\',name,'.pdf']);
%    exportgraphics(FIG1,picname1,'ContentType','vector')
   saveas(FIG1,picname1)
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
 
 POS.xlim = [-0.1,0.4];
 POS.ylim = [-1.5 , 1.5];
 
FIG2=figure(fig_num);
plot(feval(public.function.HX,C,[MINUS.mu;0],MINUS.mu),0,'ro','LineWidth', 1.5,'markerfacecolor',[1,1,1])
hold on
plot([0,0],10*[-1 1],'--','color',[0,0,0])
plot(X,Y,'-','linewidth',public.plot.S_po.linewidth,'color',public.plot.S_po.color)

% xlabel('$H(x)$','interpreter','latex')
% ylabel('$\frac{\rm{d}}{\rm{d} t}H(x)$','interpreter','latex')

xlim(POS.xlim)
ylim(POS.ylim)
set(gca,'fontsize',public.plot.gca.fontsize,...
        'fontname',public.plot.gca.fontname,...
        'linewidth',public.plot.gca.linewidth,...
        'Units','inches',...
        'position',[public.plot.gca.ws_x  public.plot.gca.ws_y ...
     public.plot.gca.width public.plot.gca.height])
    set(FIG2, 'PaperUnits', 'inches','PaperPosition', ...
    [0  0 ...
     public.plot.gca.ws_x+public.plot.gca.width  public.plot.gca.ws_y+public.plot.gca.height],...
    'PaperSize', [public.plot.gca.width public.plot.gca.height]);
    %save 
    name ='2D_BEB_case1_1LCO';
    picname2 = strcat(filepath,['\',name,'.pdf']);
    % using exportgraphycs
    % exportgraphics(FIG2,picname2,'ContentType','vector')
    % using saveas
    saveas(FIG2,picname2)
%% CASE2
%% -----------CASE 2: From admissible focus to pseudo focus --------------%
% stable sdmissible equilibrium when mu =-1;
% stable LCO when mu =1;

%------ define the system ------%
CASE = [];
CASE.num = 2;
CASE.description = 'AF2PF_NO_po';

a =-1 ;b=-1;
A = [0,1;b,a];
r= 7;
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
MINUS.description = 'addmissible unstable equilibrium and stable LCO';
MINUS.mu = -1;
fig_num =3;

% presicted by technique developed in paper
MINUS.po.starter = equi_type*[1;-60.3353];
MINUS.po.period = 3.753639396346425;
MINUS.po.Mono_p = [1;-0.164];
MINUS.po.Salt_p = [1;1.1482];
MINUS.po.description = 'Unstable LCO';
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
 
 MINUS.xlim = [-5,48];
 MINUS.ylim = [-52 , 110];
 
FIG3=figure(fig_num);
plot(feval(public.function.HX,C,[0;0],MINUS.mu),0,'ro','LineWidth', 1.5,'markerfacecolor',[0,0,0])
hold on

plot([0,0],200*[-1 1],'--','color',[0,0,0])
plot(X,Y,'linestyle',public.plot.US_po.style,'linewidth',public.plot.US_po.linewidth,'color',public.plot.US_po.color)

% xlabel('$H(x)$','interpreter','latex')
% ylabel('$\frac{\rm{d}}{\rm{d} t}H(x)$','interpreter','latex')
xlim(MINUS.xlim)
ylim(MINUS.ylim)
set(gca,'fontsize',public.plot.gca.fontsize,...
        'fontname',public.plot.gca.fontname,...
        'linewidth',public.plot.gca.linewidth,...
        'Units','inches',...
        'position',[public.plot.gca.ws_x  public.plot.gca.ws_y ...
     public.plot.gca.width public.plot.gca.height])
    set(FIG3, 'PaperUnits', 'inches','PaperPosition', ...
    [0  0 ...
     public.plot.gca.ws_x+public.plot.gca.width  public.plot.gca.ws_y+public.plot.gca.height],...
    'PaperSize', [public.plot.gca.width public.plot.gca.height]);
    %save 
       name ='2D_BEB_case2_usLCO';
    picname3 = strcat(filepath,['\',name,'.pdf']);
    % using exportgraphycs
    % exportgraphics(FIG2,picname2,'ContentType','vector')
    % using saveas
    saveas(FIG3,picname3)
% ------------------ mu =-1 ----------------------------%
equi_type =   1; 
POS = [];

POS.description = 'divergent orbit';
% presicted by technique developed in paper

POS.mu =  1;
POS.fs =fs;
POS.starter = [equi_type;1];
% post processing of the data
fig_num =4;
yt=0;
yb=0;
xr=0;
% try some plot in this case
int_2 = POS.starter;
 % outside LCO

t_2 = 6*pi;
[tout1,yout1,~,teout1,~,~]=...
    Single_DS_Impacting_Hybrid_system_integration(A,R,C,int_2,[0 t_2],fs,equi_type);
%
POS.init2 = int_2;
POS.tspan2 = t_2;


[ytop,ybot,xright,CASE,BEB2D]=trunc_plot(tout1,yout1,teout1,floor(t_2),CASE,POS,fig_num,BEB2D);

yt = max([yt,ytop]);
yb = min([yb,ybot]);
xr = max([xr,xright]);

 % ---- Plot---
 
 POS.xlim = [-5,35];
 POS.ylim = [-20 , 60];
 
FIG4=figure(fig_num);
% plot the unstable pseudo equilibrium
plot(feval(public.function.HX,C,[MINUS.mu;0],MINUS.mu),0,'ro','LineWidth', 1.5,'markerfacecolor',[1,1,1])
hold on
plot([0,0],60*[-1 1],'--','color',[0,0,0])

% xlabel('$H(x)$','interpreter','latex')
% ylabel('$\frac{\rm{d}}{\rm{d} t}H(x)$','interpreter','latex')
xlim(POS.xlim)
ylim(POS.ylim)
set(gca,'fontsize',public.plot.gca.fontsize,...
        'fontname',public.plot.gca.fontname,...
        'linewidth',public.plot.gca.linewidth,...
        'Units','inches',...
        'position',[public.plot.gca.ws_x  public.plot.gca.ws_y ...
     public.plot.gca.width public.plot.gca.height])
    set(FIG4, 'PaperUnits', 'inches','PaperPosition', ...
    [0  0 ...
     public.plot.gca.ws_x+public.plot.gca.width  public.plot.gca.ws_y+public.plot.gca.height],...
    'PaperSize', [public.plot.gca.width public.plot.gca.height]);
    %save 
    name ='2D_BEB_case2_1diverge';
    picname4 = strcat(filepath,['\',name,'.pdf']);
    % using exportgraphycs
    % exportgraphics(FIG2,picname2,'ContentType','vector')
    % using saveas
    saveas(FIG4,picname4)

 
    %% CASE3
    %------ define the system ------%
CASE = [];
CASE.num = 3;
CASE.description = 'AF2PF_po_S';

a =0.5 ;b=-1;
A = [0,1;b,a];
r= 0.5;
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
MINUS.description = 'addmissible unstable equilibrium';
MINUS.mu = -1;
fig_num =5;

% post processing of the data

yt=0;
yb=0;
xr=0;
% try some plot in this case
int_1 = zeros(2,1);
int_1(1) = 0;
int_1(2) = 0.001;
t_1 = 15*pi;
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

FIG5=figure(fig_num);
plot(feval(public.function.HX,C,[0;0],MINUS.mu),0,'ro','LineWidth', 1.5,'markerfacecolor',[1,1,1])
hold on
plot([0,0],15*[-1 1],'--','color',[0,0,0])
% xlabel('$H(x)$','interpreter','latex')
% ylabel('$\frac{\rm{d}}{\rm{d} t}H(x)$','interpreter','latex')
xlim([-1.1,xr])
set(gca,'fontsize',public.plot.gca.fontsize,...
        'fontname',public.plot.gca.fontname,...
        'linewidth',public.plot.gca.linewidth,...
        'Units','inches',...
        'position',[public.plot.gca.ws_x  public.plot.gca.ws_y ...
     public.plot.gca.width public.plot.gca.height])
    set(FIG5, 'PaperUnits', 'inches','PaperPosition', ...
    [0  0 ...
     public.plot.gca.ws_x+public.plot.gca.width  public.plot.gca.ws_y+public.plot.gca.height],...
    'PaperSize', [public.plot.gca.width public.plot.gca.height]);
    %save 
    name ='2D_BEB_case3_1diverge';
    picname5 = strcat(filepath,['\',name,'.pdf']);
    % using exportgraphycs
    % exportgraphics(FIG2,picname2,'ContentType','vector')
    % using saveas
    saveas(FIG5,picname5)
 
% ------------------ mu =-1 ----------------------------%
equi_type =   1; 
POS = [];

POS.description = 'unstable periodic orbit';
% presicted by technique developed in paper
POS.po.starter = [1;5.8733];
POS.po.period = 0.9503 * pi;
POS.po.Mono_p = [1;-2.2245];
POS.po.Salt_p = [1;1.1122];
POS.po.description = 'Stable LCO';
POS.mu =  1;
POS.fs =fs;

% post processing of the data
fig_num =6;
yt=0;
yb=0;
xr=0;
% try some plot in this case
int_2 = POS.po.starter;
 % outside LCO
int_2(2) = int_2(2) + 0.1;
t_2 = 20*pi;
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
% int_1 = POS.po.starter;
% int_1(2) = int_1(2) -3;
% t_1 = 30*pi;
% [tout1,yout1,~,teout1,~,~]=...
%     Single_DS_Impacting_Hybrid_system_integration(A,R,C,int_1,[0 t_1],fs,equi_type);
% %
% POS.init1 = int_1;
% POS.tspan1 = t_1;
% 
% [ytop,ybot,xright,CASE,BEB2D]=trunc_plot(tout1,yout1,teout1,floor(t_1),CASE,POS,fig_num,BEB2D);
% 
% yt = max([yt,ytop]);
% yb = min([yb,ybot]);
% xr = max([xr,xright]);

%
int_1 = POS.po.starter;
int_1(2) = int_1(2) -0.1;
t_1 = 100*pi;
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
 
 POS.xlim = [-1,10];
 POS.ylim = [-15 , 10];
 
FIG6=figure(fig_num);
plot(feval(public.function.HX,C,[MINUS.mu;0],MINUS.mu),0,'ro','LineWidth', 1.5,'markerfacecolor',[0,0,0])
hold on

plot([0,0],15*[-1 1],'--','color',[0,0,0])
plot(X,Y,'linestyle',public.plot.US_po.style,'linewidth',public.plot.US_po.linewidth,'color','red')

% xlabel('$H(x)$','interpreter','latex')
% ylabel('$\frac{\rm{d}}{\rm{d} t}H(x)$','interpreter','latex')

xlim(POS.xlim)
ylim(POS.ylim)
set(gca,'fontsize',public.plot.gca.fontsize,...
        'fontname',public.plot.gca.fontname,...
        'linewidth',public.plot.gca.linewidth,...
        'Units','inches',...
        'position',[public.plot.gca.ws_x  public.plot.gca.ws_y ...
     public.plot.gca.width public.plot.gca.height])
    set(FIG6, 'PaperUnits', 'inches','PaperPosition', ...
    [0  0 ...
     public.plot.gca.ws_x+public.plot.gca.width  public.plot.gca.ws_y+public.plot.gca.height],...
    'PaperSize', [public.plot.gca.width public.plot.gca.height]);
    %save 
     name ='2D_BEB_case3_2usLCO';
    picname6 = strcat(filepath,['\',name,'.pdf']);
    % using exportgraphycs
    % exportgraphics(FIG2,picname2,'ContentType','vector')
    % using saveas
    saveas(FIG6,picname6)

    
    %% CASE4
    %% -----------CASE 4: From admissible focus to pseudo focus --------------%
% stable sdmissible equilibrium when mu =-1;
% stable LCO when mu =1;

%------ define the system ------%
CASE = [];
CASE.num = 4;
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
 
FIG7=figure(fig_num);
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
        'linewidth',public.plot.gca.linewidth,...
        'Units','inches',...
        'position',[public.plot.gca.ws_x  public.plot.gca.ws_y ...
     public.plot.gca.width public.plot.gca.height])
    set(FIG7, 'PaperUnits', 'inches','PaperPosition', ...
    [0  0 ...
     public.plot.gca.ws_x+public.plot.gca.width  public.plot.gca.ws_y+public.plot.gca.height],...
    'PaperSize', [public.plot.gca.width public.plot.gca.height]);
    %save 
    name ='2D_BEB_case4_1LCO';
    picname7 = strcat(filepath,['\',name,'.pdf']);
    % using exportgraphycs
    % exportgraphics(FIG2,picname2,'ContentType','vector')
    % using saveas
    saveas(FIG7,picname7)%save
   
   
    
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
% figure(1)
% plot(yout1(:,1),yout1(:,2),'k--','linewidth',0.8)
%
POS.init = int_1;
POS.tspan = t_1;
POS.fs =fs;

[ytop,ybot,xright,CASE,BEB2D]=trunc_plot(tout1,yout1,teout1,floor(t_1),CASE,POS,fig_num,BEB2D);
yt = max([yt,ytop]);
yb = min([yb,ybot]);
xr = max([xr,xright])-equi_type+0.5;

FIG8=figure(fig_num);
plot(feval(public.function.HX,C,[1;0],POS.mu),0,'ro','LineWidth', 1.5,'markerfacecolor',[0,0,0])
hold on
plot([0,0],25*[-1 1],'--','color',[0,0,0])
% xlabel('$H(x)$','interpreter','latex')
% ylabel('$\frac{\rm{d}}{\rm{d} t}H(x)$','interpreter','latex')
xlim([-1.1,6])
ylim([-10 10])
set(gca,'fontsize',public.plot.gca.fontsize,...
        'fontname',public.plot.gca.fontname,...
        'linewidth',public.plot.gca.linewidth,...
        'Units','inches',...
        'position',[public.plot.gca.ws_x  public.plot.gca.ws_y ...
     public.plot.gca.width public.plot.gca.height])
    set(FIG8, 'PaperUnits', 'inches','PaperPosition', ...
    [0  0 ...
     public.plot.gca.ws_x+public.plot.gca.width  public.plot.gca.ws_y+public.plot.gca.height],...
    'PaperSize', [public.plot.gca.width public.plot.gca.height]);
    %save 
    name ='2D_BEB_case4_2equi';
    picname8 = strcat(filepath,['\',name,'.pdf']);
    % using exportgraphycs
    % exportgraphics(FIG2,picname2,'ContentType','vector')
    % using saveas
    saveas(FIG8,picname8)

