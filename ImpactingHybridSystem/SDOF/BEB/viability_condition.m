% viability condtion check for the algorithm

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
public.plot.gca.ws_x = 0.48*2;
public.plot.gca.ws_y = 0.3354+0.1;
public.plot.gca.width = 4.5208;
public.plot.gca.height = 3.5656;

% observing function
HX = @(C,Y,mu) C*Y - mu;
HV = @(C,A,Y) C*A*Y;
public.function.HX = HX;
public.function.HV = HV;

%
BEB2D.public =public;

%
BEB2D.public =public;

filepath = [pwd,'\figures'];
if ~exist(filepath, 'dir')
       mkdir(filepath)
end
%% -----------CASE 1: From admissible focus to pseudo focus --------------%
% stable sdmissible equilibrium when mu =-1;
% stable LCO when mu =1;
equi_type =1;
%------ define the system ------%
CASE = [];
CASE.num = 1;
CASE.description = 'AF2PF_po_S';

alpha = -0.1;
beta = 1;

lambda_1 = alpha + 1j*beta;
lambda_2 = alpha - 1j*beta;

% a =-1 ;b=-1;

a = lambda_1 + lambda_2;
b = - lambda_1 * lambda_2;


A = [0,1;b,a];

% r= 1.5;
r = 1.2*exp(-3*pi*alpha);
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

% EA = @(T) real(V1*diag(exp(diag(D1)*T))*inv(V1));
sign_V = @(Y) sign(C*A*Y);

% determine the sampling frequency 
Omega = max(abs(diag(D1)));
fs = 200*ceil(2*(Omega/2/pi));
% for different Evolution time
a =0;
b = 10;
delta = 0.01;
T = a*pi:delta:b*pi;

%
MAX = zeros(1,length(T));
F_1 = zeros(1,length(T));
V_sign = zeros(1,length(T));
LOCI= zeros(size(A,1),length(T));
vector= zeros(size(A,1),length(T));
for i=1:length(T)
    [V_sign(i),LOCI(:,i),MAX(i),vector(:,i),F_1(i)] = LCO_Det_search(T(i),R,A,C,equi_type);
end

%% make the plot
velocity = C*A*vector;FIG1 = figure;set(FIG1,'Units','inches');
hax1 = axes;

plot(T/pi,F_1,'r-','linewidth',1.2,'displayname','indication function')
hold on
%  plot( T/pi,V_sign,'-','linewidth',1.2,'displayname','Sign of Velocity')
plot([a b],[0 0],'b-','displayname','zero line')


legend
% ylim([-2 2])
xlabel('\tau /\pi')
ylabel('$p(t)$', 'interpreter','latex')
set(gca,'fontname','Times New Roman')

set(gca,'fontsize',public.plot.gca.fontsize,...
        'fontname',public.plot.gca.fontname,...
        'linewidth',public.plot.gca.linewidth,...
        'Units','inches',...
        'position',[public.plot.gca.ws_x  public.plot.gca.ws_y ...
     public.plot.gca.width public.plot.gca.height])
 pos = get(FIG1,'Position');
    set(FIG1, 'PaperUnits', 'inches','PaperPosition', ...
    [0  0 ...
    public.plot.gca.width  public.plot.gca.height],...
    'PaperSize', pos(3:4));
 %save 

index0 =sign(F_1);
index1 = abs(diff(index0))>0;
% filter the singularity case
index_s = abs(diff(F_1))/delta < (1/delta);
index1 =index1 & index_s;
index2 = [index1,0];
index3 = [0,index1];
index2=find(index2==1);
index3=find(index3==1);
% section partia
ratio = abs(F_1(index2))./(abs(F_1(index2))+abs(F_1(index3)));
T_chosen =(1-ratio).*T(index2)+ratio.*T(index3);
MAX_chosen = (1-ratio).*MAX(index2)+ratio.*MAX(index3);
F1_chosen = (1-ratio).*F_1(index2)+ratio.*F_1(index3);
Sign_chosen =(1-ratio).*V_sign(index2)+ratio.*V_sign(index3);
%          figure(1)
% plot(T_chosen/pi,Sign_chosen,'bo')
%     %
index4 =  (abs(F1_chosen)<1e-3) & (T_chosen>1e-3);
% index4 = (Sign_chosen == 1) & (abs(F1_chosen)<1e-3) %& (abs(MAX_chosen)<1e-3);
T_chosen  =T_chosen(index4);
F1_chosen =F1_chosen(index4);

if ~isempty(T_chosen)
    LCO = zeros(length(C),length(T_chosen));
    if length(T_chosen)>2
    left = T_chosen(1)-0.2*sum(abs(T_chosen(end)-T_chosen(1)));
    right =T_chosen(end)+0.2*sum(abs(T_chosen(end)-T_chosen(1)));
    else
    left = T_chosen(1)-0.1*sum(abs(T_chosen));
    right =T_chosen(end)+0.1*sum(abs(T_chosen));
    end
    for i=1:length(T_chosen)
        % get the LCO with statevariables
        [V_sign_LCO(i),~,F11,LCO(:,i),~] = LCO_Det_search(T_chosen(i),R,A,C,equi_type);
        figure(1)
        plot(T_chosen(i)/pi,F1_chosen(i),'bo','displayname',['the ',num2str(i),'th candidate'])
    end
    figure(1)
%     xlim([0 right])
    legend('location','best')
    disp([num2str(length(T_chosen)),' LCO(s) found!'])
    LCO
    disp('with period:')
    T_chosen
   
    disp('with approximation:')
     disp('with sign')
    V_sign_LCO
    
    F11
    for i=1:length(T_chosen)
        %check the floque multipliers of the foud LCO
        [Mono_p,Salt_p]=Floque_Multipliers(T_chosen(i),LCO(:,i),R,A,C)         
        if max(abs(Salt_p))>1+1e-6
            disp([' LCO',num2str(i),' is unstable!'])
        else
            disp([' LCO',num2str(i),' is stable!'])
        end
    end
else
    disp('No LCO found!')
end

keyboard

% if ~isempty(T_chosen)
%     tspan = [0 50];
%     for i =1:length(T_chosen)
%     InitCond = LCO(:,i);
%     [tout,yout,yeout0,teout,yeout,ieout]=...
%         Single_DS_Impacting_Hybrid_system_integration(A,R,C,V_sign_LCO(i)*InitCond,tspan,fs,V_sign_LCO(i)*equi_type);
%     eval(['tout',num2str(i),'=tout;']);
%     eval(['yout',num2str(i),'=yout;']);
%     eval(['yeout0',num2str(i),'=yeout0;']);
%     eval(['teout',num2str(i),'=teout;']);
%     eval(['yeout',num2str(i),'=yeout;']);
%     eval(['ieout',num2str(i),'=ieout;']);
%     eval(['delta_t',num2str(i),'=teout',num2str(i),'(2:end)-teout',num2str(i),'(1:end-1);']);
%     %Period
%     eval(['Period',num2str(i),'=delta_t',num2str(i),'/pi;']);
%     figure
%     plot(tout/pi,yout(:,C>0),'r-','linewidth',1.4)
%     title(['LCO ',num2str(i),'''s stability'])
%     set(gca,'fontname','times new roman','fontsize',12)
%     xlabel('t/s')
%     figure
%     plot(yout(:,1),yout(:,2),'k-','linewidth',1.2)
%     title(['LCO ',num2str(i),'''s phase portrait'])
%     set(gca,'fontname','times new roman','fontsize',12)
%     grid on
% %     figure
% %     plot(yout(:,1),yout(:,2),'k-','linewidth',1.2)
% %     title(['LCO ',num2str(i),'''s phase portrait'])
% %     set(gca,'fontname','times new roman','fontsize',12)
% %     grid on
% %     figure
% %     plot(yout(:,2),yout(:,3),'k-','linewidth',1.2)
% %     title(['LCO ',num2str(i),'''s phase portrait'])
% %     set(gca,'fontname','times new roman','fontsize',12)
% %     grid on
%     end
% end
% %% 
InitCond = LCO(:,1);
    [tout,yout,yeout0,teout,yeout,ieout]=...
        Single_DS_Impacting_Hybrid_system_integration(A,R,C,V_sign_LCO(1)*InitCond,[0 T_chosen(1)],fs,100*V_sign_LCO(1)*equi_type);
 figure
    plot(yout(:,1),yout(:,2),'k-','linewidth',1.6)
    hold on
    plot(-[1,1],200*[-1 1],'--','color',[0,0,0])
    ylim([-2. 2.])
    xlim([-2 2])
InitCond = LCO(:,2);
    [tout,yout,yeout0,teout,yeout,ieout]=...
        Single_DS_Impacting_Hybrid_system_integration(A,R,C,V_sign_LCO(2)*InitCond,[0 T_chosen(2)],fs,100*V_sign_LCO(2)*equi_type);
 figure
    plot(yout(:,1),yout(:,2),'k-','linewidth',1.6)
    hold on
    plot([1,1],200*[-1 1],'--','color',[0,0,0])
    ylim([-2.5 2.5])
%%
InitCond = LCO(:,3);
    [tout,yout,yeout0,teout,yeout,ieout]=...
        Single_DS_Impacting_Hybrid_system_integration(A,R,C,V_sign_LCO(3)*InitCond,[0 T_chosen(3)],fs,100*V_sign_LCO(3)*equi_type);
 figure
    plot(yout(:,1),yout(:,2),'k-','linewidth',1.6)
    hold on
    plot(-[1,1],200*[-1 1],'--','color',[0,0,0])
    ylim([-6 8])