% Detecting the existence of LCO when boundary equilibrium transferring to
% virtual Equilibrium for general n dimension impacting hybrid system 
% limitations: single impact; single discontinuity surface

% KEY INPUTs:
% A : the linear region system's matrix 
% R : the reset map matrix, R(x) = R*x;
% C : the selecting vector to define the discontinuity surface C*x=0
% 

% TRY with a 3D problem
clc
clear
close all

rootpath = erase(mfilename('fullpath'),mfilename());
filepath = [rootpath,'\figures'];
if ~exist(filepath, 'dir')
       mkdir(filepath)
end

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
% if equi_type == 'admissible'
%         pre = -1;
%  elseif equi_type == 'pseudo'
%         pre = 1;
%  end
equi_type = 1;

% 
% define the linear part
alpha = -0.1;
beta  = 0.2;
lambda = -0.5;
A_1=[alpha beta 0; -beta alpha 0;0 0 lambda];
diag(eig(A_1))
[V,D]=eig(A_1)
theta = pi/6;
P = [cos(theta),0,-sin(theta);0,1,0;sin(theta),0,cos(theta)];
A = inv(P)*A_1*P;

% define the reset map
a11=A(1,1);a12=A(1,2);a13=A(1,3);
N = sqrt(a12^2+a13^2);
V_l3 = [0;a13;-a12]/N/N;
P_l3 = [0;a12; a13]/N/N;

pr= 0.8; %coef of restitution
pz = 4.6;

MW =  -(1+pr)*P_l3*[a11 a12 a13]+pz*V_l3*[a11 a12 a13];
B = @(P_,V_) -(1+pr)*P_+pz*V_;

% define initial condition: grid on the positive set
C = [1 0 0];
R = eye(length(C))+MW;


%% define the system in the form of Lenard form
% A = [-0.7 1 0;-0.15 0 1;-0.025 0 0]
% B =[0;2.5;0.625]
% C = [1 0 0];
% R = eye(3) - B*C*A
[V1,D1]=eig(A);
EA = @(T) real(V1*diag(exp(diag(D1)*T))*inv(V1));
sign_V = @(Y) sign(C*A*Y);
% determine the sampling frequency 
Omega = max(abs(diag(D1)));
fs = 10*ceil(2*(Omega/2/pi));
% for different Evolution time
a =0;
b = 3;
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
ylabel('value')
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
    name ='3D_case1_alg';
    picname1 = strcat(filepath,['\',name,'.pdf']);
    
%%
B0 =B(P_l3,V_l3)

[A_c,B_c,PZ] = trans2companion_form(A,-B0,B0)


[A_c,B_c,PZ] = trans2companion_form(A,-B0,P_l3)

eig(A_c)
B_c = B((PZ\P_l3)/norm(PZ\P_l3),(PZ\V_l3)/norm(PZ\V_l3))
% stability of PE
A_s_c = (eye(3)-(B_c*C*A_c)/(C*A_c*B_c))*A_c;
asc =eig(A_s_c)

A_s = (eye(3)-(B0*C*A)/(C*A*B0))*A;
as=eig(A_s)
pec = pinv([C;C*A_c;A_s_c])*[C';0;0]
pe = pinv([C;C*A;A_s])*[C';0;0]

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
MAX_chosen  =MAX_chosen(index4);
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
        [~,~,~,LCO(:,i),~] = LCO_Det_search(T_chosen(i),R,A,C,equi_type);
        figure(1)
        plot(T_chosen(i)/pi,F1_chosen(i),'bo','displayname',['the ',num2str(i),'th candidate'])
    end
    figure(1)
    xlim([0 right]/pi)
    legend('location','best')
    disp([num2str(length(T_chosen)),' LCO(s) found!'])
    LCO
    disp('with period:')
    T_chosen/pi
    disp('with approximation:')
    
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
%
if ~isempty(T_chosen)
    tspan = [0 50];
    for i =1:length(T_chosen)
    InitCond = LCO(:,i);
    [tout,yout,yeout0,teout,yeout,ieout]=...
        Single_DS_Impacting_Hybrid_system_integration(A,R,C,InitCond,tspan,fs,equi_type);
    eval(['tout',num2str(i),'=tout;']);
    eval(['yout',num2str(i),'=yout;']);
    eval(['yeout0',num2str(i),'=yeout0;']);
    eval(['teout',num2str(i),'=teout;']);
    eval(['yeout',num2str(i),'=yeout;']);
    eval(['ieout',num2str(i),'=ieout;']);
    eval(['delta_t',num2str(i),'=teout',num2str(i),'(2:end)-teout',num2str(i),'(1:end-1);']);
    %Period
    eval(['Period',num2str(i),'=delta_t',num2str(i),'/pi;']);
    figure
    plot(tout/pi,yout(:,C>0),'r-','linewidth',1.4)
    title(['LCO ',num2str(i),'''s stability'])
    set(gca,'fontname','times new roman','fontsize',12)
    xlabel('t/s')
    figure
    plot(yout(:,1),yout(:,2),'k-','linewidth',1.2)
    title(['LCO ',num2str(i),'''s phase portrait'])
    set(gca,'fontname','times new roman','fontsize',12)
    grid on
    end
end
