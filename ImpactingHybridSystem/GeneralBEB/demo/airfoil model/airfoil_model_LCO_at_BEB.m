addpath ../public_scripts
% 初始化
clc
close all
clear

%% define the system
Stiffness_ratio=1;
% start computation
gap=0.01;
heave0=0;
alpha0=0.06;
beta0=1e-4/Stiffness_ratio;
damp0=1;
%
Stiff_k1=1;
Stiff_k2=300*sqrt(Stiffness_ratio);

method='Pade';
FLUTEER_MATRIXBASE=Initialization_backlash(gap,heave0,alpha0,beta0,damp0,Stiff_k1,Stiff_k2);
% 定义间隙大小
gap=FLUTEER_MATRIXBASE.TimeInt.gap;
MHH=FLUTEER_MATRIXBASE.MHH;
fs=FLUTEER_MATRIXBASE.TimeInt.fs;
Init=FLUTEER_MATRIXBASE.Init;
%
V_care=0.64833*30;
res=0.72;
[A_1,A_2,M_new,K_new,C_new,vect_preload]=Get_int_Matrix(FLUTEER_MATRIXBASE,method,V_care);
[V,D]=eig(A_1);
Omega=diag(D)
%% define the reset map & get the matrix of P
BEquilibrium=Static_equilibrium_output(gap,A_1,A_2,vect_preload);
% InitCond(3)=-0.03;
%   InitCond(7)=0;
 [A,M,C,N,B]=Get_Impact_BEB_Matrix(FLUTEER_MATRIXBASE,method,V_care,BEquilibrium,res);
 % composed map Jacoi matrix and its eigenvalues
MW=zeros(8);
MW(6,6)=B(6);
MW(4,6)=B(4);
MW(5,6)=B(5);
P=eye(8)+MW;
% different time length for the decaying process
T=100;
dt=1/fs;
t=0:dt:T;

%% 时域方法： 一阶数值积分
num=44;
tspan=[0 T];
n_init=5;
Balance_state1=-A_1\vect_preload(:,1)
Balance_state2=-A_2\vect_preload(:,1);
scaler = abs(Balance_state1(3)- -gap);
LCO1=[-0.130979630614067;-0.0157740175410459;1;8.56000986141258;-5.66866590417186;69.5166411243114;1.56431285757849;-0.255799411131381];
viation = LCO1*scaler;
% InitCond = [heave0;alpha0;gap*(-1+(abs(n_init)-1)/(0.25*num-0.5));0;0;0*20*gap*sign(n_init);0;0];
InitCond = [-0.00463900035393225,0.0616328862516022,-0.00991614723523688,-0.000235257200054251,-0.000175628954363620,0.00121299515977076,0.0623356166075342,-3.12652396655737e-06]';

%  InitCond(6,1) = -20*gap*sign(InitCond(3,1));
%  InitCond(2,1) = 0;
%  InitCond(6,1) = 0.0;

correction=Balance_state1;
% correction(3)=-0.01;
% InitCond=Static_equilibrium_output(gap,A_1,A_2,vect_preload)
% %
% vect_preload(:,1)=A_1*(correction-Balance_state1)
% A_1(2:8,2:8)\vect_preload(2:8,1)

%
% InitCond=viation+Balance_state1
% InitCond(3)=-0.01;
% InitCond(6)=0.1;
%  InitCond(7)=0;
%   InitCond(8)=0;
% backlash_airfoil_restitution(A_1,A_2,InitCond,tspan);
% [tout,yout,teout,yeout,ieout,amplitude]=backlash_airfoil_K(A_1,A_2,InitCond,vect_preload,tspan,gap,fs);


% InitCond=impact_map(MHH,res,InitCond);
[tout,yout,yeout0,teout,yeout,ieout,amplitude]=backlash_airfoil_restitution(res,A_1,MHH,InitCond,vect_preload,tspan,gap,fs,1);
% [tout,yout,yeout0,teout,yeout,ieout,amplitude]=backlash_airfoil_restitution_H2(res,A_1,MHH,InitCond,vect_preload,tspan,gap,fs);
% [tout,yout,yeout0,teout,yeout,ieout]=Composed_evolution_airfoil(A_1,P,1.01*LCO1,tspan,fs);
delta_t=teout(2:end)-teout(1:end-1);
%Period
Period=teout(2:end)-teout(1:end-1);
Period=Period/pi;

%% Theoritical balance state
Preload =FLUTEER_MATRIXBASE.Preload;
K1      =FLUTEER_MATRIXBASE.Khh;
K2      =FLUTEER_MATRIXBASE.KHH;
MAC     =FLUTEER_MATRIXBASE.MAC;
b       =0.5*MAC;
Q_d     =FLUTEER_MATRIXBASE.Q_d;
A0=FLUTEER_MATRIXBASE.Pade.A0;
Kernel1=K1-Q_d*(V_care/b)^2*A0;
Kernel2=K2-Q_d*(V_care/b)^2*A0;
% state=vect_preload\A_2;
% Balance_state1=Kernel1\Preload;
% Balance_state2=Kernel2\Preload;


%% Time series Plot
 %% beta motion
%         find v==0 points with set tol ----poincare section method I

% heave motion
figure(1)
plot(tout,yout(:,1),'r-')
hold on
if ~isempty(yeout)
plot(teout,yeout(:,1),'ko');
end
% if ~isempty(yeout)
%     hold on
%     plot(tout(IndMin),yout(IndMin,3),'ko')
%     plot(tout(IndMax),yout(IndMax,3),'ko')
% end
legend('falp')
xlabel('time');
ylabel('dis');
title('heave motion with backlash');
set(gca,'fontname','Times New Roman','fontsize',12)
%pitch motion
figure(2)
plot(tout,yout(:,2),'r-')
hold on
 if ~isempty(yeout)
plot(teout,yeout(:,2),'ko');
 end
% if ~isempty(yeout)
%     hold on
%     plot(tout(IndMin),yout(IndMin,3),'ko')
%     plot(tout(IndMax),yout(IndMax,3),'ko')
% end
legend('falp')
xlabel('time');
ylabel('relative angle/rad');
title('pitch motion with backlash');
set(gca,'fontname','Times New Roman','fontsize',12)
%
figure(3)
plot(tout,yout(:,3),'r-')
hold on
 if ~isempty(yeout)
plot(teout,yeout(:,3),'ko');
 end
% if ~isempty(yeout)
%     hold on
%     plot(tout(IndMin),yout(IndMin,3),'ko')
%     plot(tout(IndMax),yout(IndMax,3),'ko')
% end
legend('falp')
xlabel('time');
ylabel('relative angle/rad');
title('flap motion with backlash');
set(gca,'fontname','Times New Roman','fontsize',12)
%
figure(4)
plot(tout,yout(:,4),'k-')
hold on
data=yout(:,4)';
% IndMin=find(diff(sign(diff(data)))>0)+1;   %locally minumum index
% IndMax=find(diff(sign(diff(data)))<0)+1;   %locally maximum index
% if ~isempty(yeout)
%     hold on
%     plot(tout(IndMin),yout(IndMin,6),'ro')
%     plot(tout(IndMax),yout(IndMax,6),'ro')
% end
 if ~isempty(yeout)
plot(teout,yeout(:,4),'ro');
 end
legend('heave')
xlabel('time');
ylabel('relative angle velocity/rad');
title('heave motion with backlash');
set(gca,'fontname','Times New Roman','fontsize',12)
% 
figure(5)
plot(tout,yout(:,5),'k-')
hold on

% IndMin=find(diff(sign(diff(data)))>0)+1;   %locally minumum index
% IndMax=find(diff(sign(diff(data)))<0)+1;   %locally maximum index
% if ~isempty(yeout)
%     hold on
%     plot(tout(IndMin),yout(IndMin,6),'ro')
%     plot(tout(IndMax),yout(IndMax,6),'ro')
% end 
if ~isempty(yeout)
plot(teout,yeout(:,5),'ro');
end
legend('falp')
xlabel('time');
ylabel('pitch velocity/rad');
title('pitch motion with backlash');
set(gca,'fontname','Times New Roman','fontsize',12)
% plot the flap velocity 
figure(6)
plot(tout,yout(:,6),'k-')
hold on
data=yout(:,6)';
% IndMin=find(diff(sign(diff(data)))>0)+1;   %locally minumum index
% IndMax=find(diff(sign(diff(data)))<0)+1;   %locally maximum index
% if ~isempty(yeout)
%     hold on
%     plot(tout(IndMin),yout(IndMin,6),'ro')
%     plot(tout(IndMax),yout(IndMax,6),'ro')
% end
 if ~isempty(yeout)
plot(teout,yeout(:,6),'ro');
 end
legend('falp')
xlabel('time');
ylabel('relative angle velocity/rad');
title('flap motion with backlash');
set(gca,'fontname','Times New Roman','fontsize',12)

%
figure(7)
plot(tout,yout(:,7),'k-')
% if ~isempty(yeout)
%     hold on
%     plot(tout(IndMin),yout(IndMin,7),'ro')
%     plot(tout(IndMax),yout(IndMax,7),'ro')
% end
hold on
 if ~isempty(yeout)
plot(teout,yeout(:,7),'ro');
 end
legend('aero motion 1')
xlabel('time');
ylabel('relative angle/rad');
title('flap motion with backlash');
set(gca,'fontname','Times New Roman','fontsize',12)
%
%
figure(8)
plot(tout,yout(:,8),'k-')
% if ~isempty(yeout)
%     hold on
%     plot(tout(IndMin),yout(IndMin,7),'ro')
%     plot(tout(IndMax),yout(IndMax,7),'ro')
% end
hold on
 if ~isempty(yeout)
plot(teout,yeout(:,8),'ro');
 end
legend('aero motion 2')
xlabel('time');
ylabel('relative angle/rad');
title('flap motion with backlash');
set(gca,'fontname','Times New Roman','fontsize',12)
%
figure
subplot(311)
plot(tout,yout(:,1),'r-','linewidth',1.2,'displayname','heave motion')
legend
% xlim([99 100])
subplot(312)
plot(tout,yout(:,2),'k-','linewidth',1.2,'displayname','pitch motion')
legend
% xlim([99 100])
subplot(313)
plot(tout,yout(:,3),'g-','linewidth',1.2,'displayname','flap motion')
legend

%
figure
plot(yout(:,3),yout(:,6))
hold on
plot(InitCond(3),InitCond(6),'r*')
% phase portrait of Beta motion
figure
plot(yout(end-5*fs:end,3),yout(end-5*fs:end,6),'r.')
xlabel('\beta');
ylabel('$\dot{\beta}$','Interpreter','latex');
title('flap motion phase portrait');
% % xlim([-0.012 0.012])
% % ylim([-0.12 0.12])
% set(gca,'fontname','Times New Roman','fontsize',12)
% hold off
% 
% % phase portrait of pitch motion
% figure(10)
% plot(yout(end-5*fs:end,2),yout(end-5*fs:end,5),'r.')
% xlabel('\beta');
% ylabel('$\dot{\beta}$','Interpreter','latex');
% title('pitch motion phase portrait');
% % xlim([-0.012 0.012])
% % ylim([-0.12 0.12])
% set(gca,'fontname','Times New Roman','fontsize',12)
% hold off
%
% figure
% plot(yout(1:0.4*fs,3),yout(1:0.4*fs,6),'r.')
% xlabel('\beta');
% ylabel('$\dot{\beta}$','Interpreter','latex');
% title('flap motion phase portrait');
% xlim([-0.012 0.012])
% ylim([-0.12 0.12])
% set(gca,'fontname','Times New Roman','fontsize',12)
% hold off

%% Spectrum analysis
y_trunct=yout(end-20*fs:end,7);
len=length(y_trunct);

f1=(0:len/2)*fs/len;
Sx=fft(y_trunct);
P1=2*abs(Sx(1:len/2+1))/len;%transform to one-side spectrum
figure
plot(f1,P1,'k-','Displayname','frequency','linewidth',1.4)
xlim([1 80])
xlabel('\itFrequency/Hz')
ylabel('\itAmplitude')
title('Spectrum')
% legend1=legend('show','location','best');
% set(legend1,'Position',[0.25 0.8 0.12 0.05]);
set(gca,'fontname','times new roman','fontangle','italic','fontsize',14)
% % 
% %% Time-Frequency analysis
% Stfft_application
% plot the flap velocity 

select=[1,2,3,4,5,6];
tangent=A_1(:,select)*yout(:,select)'+1*vect_preload(:,1);
figure(12);
plot(tout,tangent(6,:))
title('the acceleration of flap')
figure(13)
plot(tout,tangent(6,:)/max(tangent(6,:)))
hold on
plot(tout,yout(:,3)+0.01)
figure(14)
plot(tout(end-10*fs:end),(yout(end-10*fs:end,7)-0.5*(max(yout(end-10*fs:end,7))+min(yout(end-10*fs:end,7))))/(0.5*(max(yout(end-10*fs:end,7))-min(yout(end-10*fs:end,7)))),'k-')
hold on
plot(tout(end-10*fs:end),(yout(end-10*fs:end,8)-0.5*(max(yout(end-10*fs:end,8))+min(yout(end-10*fs:end,8))))/(0.5*(max(yout(end-10*fs:end,8))-min(yout(end-10*fs:end,8)))),'r-')
plot(tout(end-10*fs:end),(yout(end-10*fs:end,3)+0.01)/max(yout(end-10*fs:end,3)+0.01))
%% observe the Kinetic energy change 
yeout(end,4:6)*M_new*yeout(end,4:6)'
yeout0(end,4:6)*M_new*yeout0(end,4:6)'
for kk=1:length(tout)
    
Kinetic(kk)=yout(kk,4:6)*M_new*yout(kk,4:6)';
Potential(kk)=yout(kk,1:3)*K_new*yout(kk,1:3)';
end
figure
subplot(311)
plot(tout,Kinetic,'k-','linewidth',1.2)
% xlim([99 100])
subplot(312)
plot(tout,Potential,'k-','linewidth',1.2)
% xlim([99 100])
subplot(313)
plot(tout,Kinetic+Potential,'k-','linewidth',1.2)
% xlim([99 100])
function y0=impact_map(MHH,res,y0)
delta_v=MHH(1:2,1:2)\[MHH(1,3);MHH(2,3)]*(1+res)*y0(6);
y0(6)=-res.*y0(6);
y0(4)=y0(4)+delta_v(1);
y0(5)=y0(5)+delta_v(2);
end



function [tout,yout,yeout0,teout,yeout,ieout]=Composed_evolution_airfoil(A_1,P,InitCond,tspan,fs)
%*************************#####************************* #
% checked and modified by Peter Tang /11th,3,2021       *#
%#************************#####************************* #
% gap:  backlash amount
% define integration timestep and timespan
tstart = tspan(1);
tfinal = tspan(2);

% define initial condition to start
y0 = InitCond;
% % check if start at the boundary
% if abs(y0(3)^2/gap^2-1)<1e-3 && y0(3)*y0(6)>0
%     y0=impact_map(MHH,res,y0);
% end
BeStuck=false;
% define the matrix of this linear differential system
A_matrix=A_1;
%
f_ls= @(t,y) A_matrix*y;
% define options
refine = 1;
options = odeset('Events',@(t,y) DetectingContact_events(t,y),...
    'RelTol',1e-12,'AbsTol',1e-12,'Refine',refine);
opsl_stuck=odeset(options,'Events',@(t,y) ...
    unstickEvent(t,y,A_matrix,vect_preload,mu));

% start calculation
tout = tstart;
yout = y0.';
teout = [];
yeout = [];
yeout0 = [];
ieout = [];

% initialization for different contact mode
Event_flag=[];

while tout(end)<tfinal
    % define out put time serias
    try
        timespan=[tstart:1/fs:tfinal];
    catch ME
        keyboard
    end
    %     timespan=[tstart tfinal];
    
    % exit case I
    if length(timespan)<2
        % last time step
        break;
    end
    
    if BeStuck==false
        % Solve until the first terminal event: hit the boundary.
        [t,y,te,ye,ie] = ode45(@(t,y) f_ls(t,y),...
            timespan,y0,options);
        %         options = odeset(options,'InitialStep',t(end)-t(end-1));
        options = odeset(options,'InitialStep',1e-6,'MaxStep',1e-4);
        %         disp('hit the boundary one time')
        Event_flag='hit boundary';
    end
    % Accumulate output.  This could be passed out as output arguments.
    nt    = length(t);
%     if t(end)>2
%         keyboard
%     end
 
    % output history
    tout  = [tout; t(2:nt)];
    yout  = [yout; real(y(2:nt,:))];
    
    % output events
    teout = [teout; te];
    yeout = [yeout; ye];
    ieout = [ieout; ie];
    if size(y0,1)>1
    yeout0 = [yeout0; y0'];
    else
        yeout0 = [yeout0; y0];
    end
    %evaluate new initial condition for next integration step if possible
    if isempty(ie)
        % means that no event detected and the integration finished till
        % the time destination, so just end the intgration and return
        % result to the mainfunction
        break
    elseif BeStuck == true
        % if previous integration is sliding, then it shows the sliding
        % stoped, continue integration
        y0    = y(nt,:);
        tstart = t(nt);
        BeStuck = false; % stuck over
        continue
    end
%    dsiP(impact_map)
    %------------------- else applay impact map--------------------------%
    t0=t(end);
    x=y(end,:)';
    y0=P*x;
    tstart=t0+0.1/fs;
    tout=[tout;tstart];
    yout=[yout;y0'];
end
%------------------------------------------------------------------------%
function [value,isterminal,direction] = DetectingContact_events(t,y)
% Locate the time when flap hit through boundary in both directions
% and stop integration when isterminal ==1
% Events: switch boundary; zero heave velocity; zero pitch velocity; zero
% flap velocity
value = [(y(3)-1);abs(y(6))];                %  detect switching point
isterminal = [1;0];                              %  stop the integration
direction  = [-1;0];                             %  mark with event funct-
                  
end
end
%%
% define the function to find the corresponding matrixes
 function [A,M,C_t,D,B]=Get_Impact_BEB_Matrix(FLUTEER_MATRIXBASE,method,V,q_0,res)
%*************************#####************************* #
% checked and modified by Peter Tang /30th,4,2021       *#
%#************************#####************************* #
% Return the state matrix for appointed method to accomplish time integration 
% process

K1      =FLUTEER_MATRIXBASE.Khh;
K2      =FLUTEER_MATRIXBASE.KHH;
MHH     =FLUTEER_MATRIXBASE.MHH;
CHH1     =FLUTEER_MATRIXBASE.CHH1;
CHH2     =FLUTEER_MATRIXBASE.CHH2;
MAC     =FLUTEER_MATRIXBASE.MAC;
NMODE   =FLUTEER_MATRIXBASE.NMODE;
b       =0.5*MAC;
Q_d     =FLUTEER_MATRIXBASE.Q_d;
Preload =FLUTEER_MATRIXBASE.Preload;


switch method
    
    case 'Pade'
%% 输出方法三：Method III: From Ref. Edwards, John W. Unsteady aerodynamic modeling for
%%						 arbitrary motions,p.372
ASE_MATRIXBASE= FLUTEER_MATRIXBASE.Pade;
A0=ASE_MATRIXBASE.A0;
A1=ASE_MATRIXBASE.A1;
A2=ASE_MATRIXBASE.A2;

S_c1=ASE_MATRIXBASE.S_c1;
S_c2=ASE_MATRIXBASE.S_c2;
R_c=ASE_MATRIXBASE.R_c;

N=size(MHH,1);
M_=MHH-Q_d*A2;

K_1=K1-Q_d*(V/b)^2*A0;  % region 1
K_2=K2-Q_d*(V/b)^2*A0;  % region 2
dot_K_V=-Q_d*2*V/b^2*A0;

B_1=CHH1-Q_d*(V/b)*A1;
B_2=CHH2-Q_d*(V/b)*A1;
D_=Q_d*(V/b)*R_c*[0.006825*(V/b)^2,0.10805*V/b];
F_p=[0 1;-0.01365*(V/b)^2 -0.3455*(V/b)];
E_1=(V/b)*[zeros(1,N);S_c1];
E_2=[zeros(1,N);S_c2];
%
dot_B_V=-Q_d/b*A1;
dot_D_V=Q_d*R_c*[3/b*0.006825*(V/b)^2,2/b*0.10805*V/b];
dot_Fp_V=[0 0;-0.01365*2*V/b^2 -0.3455/b];
dot_E1_V=[zeros(1,N);S_c1]/b;
dot_E2_V=zeros(size(E_2));
% region 1
[A]=[zeros(N),eye(N),zeros(N,2);
    -M_\K_1,-M_\B_1,M_\D_;
    E_1,E_2,F_p];

% e2=0.5*w_b_^2*r2_b_*delta*sign(beta);  % 数值计算确定beta时可用
dot_A_V=[zeros(N),0*eye(N),zeros(N,2);
    -M_\dot_K_V,-M_\dot_B_V,M_\dot_D_V;
    dot_E1_V,dot_E2_V,dot_Fp_V];
%
[M_new]=M_;
[K_new]=K_1;
[C_new]=B_1;
%%
% 
% $$e^{\pi i} + 1 = 0$$
% 
end
%%
 C_1= (MHH(2,2)*MHH(1,3)-MHH(1,2)*MHH(2,3))/(MHH(1,1)*MHH(2,2)-MHH(1,2)^2);
 C_2= (MHH(1,1)*MHH(2,3)-MHH(1,2)*MHH(1,3))/(MHH(1,1)*MHH(2,2)-MHH(1,2)^2);

M=dot_A_V*q_0;
C_t=[0;0;1;zeros(5,1)]';
D=0;
B=(1+res)*[0;0;0;C_1;C_2;-1;0;0];
end