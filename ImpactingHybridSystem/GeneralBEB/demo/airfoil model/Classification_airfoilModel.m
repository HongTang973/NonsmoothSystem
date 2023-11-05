% Using the routine to check the existence for the airfoil model Peter
% 17/02/2022


% show existence of the LCO define the system
clc
clear
close all
% if equi_type == 'admissible'
%         pre = -1;
%  elseif equi_type == 'pseudo'
%         pre = 1;
%  end
equi_type = 1;
%%  define the linear part & get the matrix A
Stiffness_ratio=1;
% start computation
gap=0.01;
heave0=0;
alpha0=0.06;
beta0=1e-4/Stiffness_ratio;
damp0=1;
%
Stiff_k1=0;
Stiff_k2=300*sqrt(Stiffness_ratio);

method='Pade';
FLUTEER_MATRIXBASE=Initialization_backlash(gap,heave0,alpha0,beta0,damp0,Stiff_k1,Stiff_k2);
% 定义间隙大小

gap=FLUTEER_MATRIXBASE.TimeInt.gap;
MHH=FLUTEER_MATRIXBASE.MHH;
fs=FLUTEER_MATRIXBASE.TimeInt.fs;
V_care=0.64833*30;
% V_care=2*30;
[A_1,A_2,M_new,K_new,C_new,vect_preload]=Get_int_Matrix(FLUTEER_MATRIXBASE,method,V_care);

res=0.72;

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
R=eye(8)+MW;
R'*R -eye(8)
eig(R)
[V1,D1]=eig(A_1);
Omega = max(abs(diag(D1)));
fs = 10*ceil(2*(Omega/2/pi));
EA = @(T) real(V1*diag(exp(diag(D1)*T))*inv(V1));
sign_V = @(Y) sign(C*A_1*Y);
%
A_s = (eye(8)-(B*C*A)/(C*A*B))*A;
pe = real(pinv([C;C*A;A_s])*[1;0;0;0;0;0;0;0;0;0])
latex(sym(pe))
eig(A_s)/100;
(A_s-A)/10000;
[V2, J]=jordan(A_1)
KK=V2*inv(R)*inv(V2)
AA=A*A
%% do line search of the time T
% for different Evolution time
a =0.0;
b = 3;
delta = 0.00005;
T = a*pi:10*delta:b*pi;
%  T = 0.0:0.01:0.1;
MAX = zeros(1,length(T));
F_1 = zeros(1,length(T));
V_sign = zeros(1,length(T));
LOCI= zeros(size(A_1,1),length(T));
vector= zeros(size(A_1,1),length(T));
for i=1:length(T)
    [V_sign(i),LOCI(:,i),MAX(i),vector(:,i),F_1(i)] = LCO_detecting_line_search(T(i),R,EA,sign_V,C,equi_type);
end

%% make the plot
figure(1)
plot( T/pi,F_1,'r-','linewidth',1.2,'displayname','F1')
hold on
% plot( T/pi,V_sign,'-','linewidth',1.2,'displayname','Sign of Velocity')
plot([a b],[0 0],'b-','displayname','zero line')
% plot( T/pi,MAX,'k-','linewidth',1.0,'displayname','MAX-1')

% plot( T/pi,LOCI(1,:),'-','linewidth',1.4,'displayname','EIG1-1') 
% plot(T/pi,LOCI(2,:),'-','linewidth',1.4,'displayname','EIG1-2')
% plot(T/pi,LOCI(3,:),'-','linewidth',1.4,'displayname','EIG1-3')
% plot( T/pi,LOCI(4,:),'-','linewidth',1.4,'displayname','EIG1-1') 
% plot(T/pi,LOCI(5,:),'-','linewidth',1.4,'displayname','EIG1-2')
% plot(T/pi,LOCI(6,:),'-','linewidth',1.4,'displayname','EIG1-3')
% plot( T/pi,LOCI(7,:),'-','linewidth',1.4,'displayname','EIG1-1') 
% plot(T/pi,LOCI(8,:),'-','linewidth',1.4,'displayname','EIG1-2')

legend
ylim([-2 2])
xlabel('\tau /\pi')
ylabel('value')
set(gca,'fontname','Times New Roman')

%% post processing

% there is a/some cross points
index0 =sign(F_1);
index1 = abs(index0(2:end)-index0(1:end-1))>0;
index2 = [index1,0];
index3 = [0,index1];
index2=find(index2==1);
index3=find(index3==1);
% section partia
ratio = abs(F_1(index2))./(abs(F_1(index2))+abs(F_1(index3)));
T_chosen =(1-ratio).*T(index2)+ratio.*T(index3);
MAX_chosen = (1-ratio).*MAX(index2)+ratio.*MAX(index3);
Sign_chosen =(1-ratio).*V_sign(index2)+ratio.*V_sign(index3);
%          figure(1)
% plot(T_chosen/pi,Sign_chosen,'bo')
%     %
% index4 = (abs(MAX_chosen)<1e-2);
index4 = (Sign_chosen == 1) & (abs(MAX_chosen)<1e-3);
T_chosen  =T_chosen(index4);
MAX_chosen =MAX_chosen(index4);
%

if ~isempty(T_chosen)
    LCO = zeros(length(C),length(T_chosen));
    left = T_chosen(1)-0.1*sum(abs(T_chosen));
    right =T_chosen(end)+0.1*sum(abs(T_chosen));
    for i=1:length(T_chosen)
        % get the LCO with statevariables
        [~,~,~,LCO(:,i),~] = LCO_detecting_line_search(T_chosen(i),R,EA,sign_V,C,equi_type);
        figure(1)
        plot(T_chosen(i)/pi,MAX_chosen(i),'bo','displayname',['the ',num2str(i),'th candidate'])
    end
    figure(1)
    xlim([left right]/pi)
    disp([num2str(length(T_chosen)),' LCO(s) found!'])
    LCO
    disp('with period:')
    T_chosen/pi
    for i=1:length(T_chosen)
        %check the floque multipliers of the foud LCO
        [Mono_p,Salt_p]=Floque_Multipliers(T_chosen(i),LCO(:,i),R,A_1,C)
        if max(abs(Salt_p))>1+1e-4
            disp([' LCO',num2str(i),' is unstable!'])
        else
            disp([' LCO',num2str(i),' is stable!'])
        end
    end
else
    disp('No LCO found!')
end
%
keyboard
%
if ~isempty(T_chosen)
    tspan = [0 20];
    for i =1:length(T_chosen)
    InitCond = LCO(:,i);
%     InitCond(3)=1;
%     InitCond(6)=0.00001;
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
    plot(tout,yout(:,C>0),'r-','linewidth',1.4)
    title(['LCO ',num2str(i),'''s stability'])
    set(gca,'fontname','times new roman','fontsize',12)
    xlabel('t/s')
    figure
    plot(yout(:,3),yout(:,6),'k-','linewidth',1.2)
    title(['LCO ',num2str(i),'''s phase portrait'])
    set(gca,'fontname','times new roman','fontsize',12)
    grid on
    end
end

Period = teout(2:end)-teout(1:end-1);
% obeserve the LCO in phase portrait
myout = yout';
figure 
subplot(221)
plot(myout(1,:),myout(4,:))
subplot(222)
plot(myout(2,:),myout(5,:))
subplot(223)
plot(myout(3,:),myout(6,:))
subplot(224)
plot(myout(7,:),myout(8,:))
% Observe the LCO in modal coordinates
myout = inv(V1)*yout';
% 
figure 
subplot(221)
plot(myout(1,:),myout(4,:))
subplot(222)
plot(myout(2,:),myout(5,:))
subplot(223)
plot(myout(3,:),myout(6,:))
subplot(224)
plot(myout(7,:),myout(8,:))



% define the function to find the corresponding matrixes
function [A,M,C_t,D,B]=Get_Impact_BEB_Matrix(FLUTEER_MATRIXBASE,method,V,q_0,res)
%*************************#####************************* #
% checked and modified by Peter Tang /30th,4,2021       *#
%#************************#####************************* #
% Return the state matrix for appointed method to accomplish time
% integration process

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