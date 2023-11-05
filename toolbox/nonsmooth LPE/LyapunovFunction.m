% track the Lyapunov function during the evolution 
addpath ../public_scripts
% 初始化
clc
close all
clear

% 导入模型数据 得到基本矩阵

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
Khh=FLUTEER_MATRIXBASE.Khh;
C_1= (MHH(2,2)*MHH(1,3)-MHH(1,2)*MHH(2,3))/(MHH(1,1)*MHH(2,2)-MHH(1,2)^2);
C_2= (MHH(1,1)*MHH(2,3)-MHH(1,2)*MHH(1,3))/(MHH(1,1)*MHH(2,2)-MHH(1,2)^2);
fs=FLUTEER_MATRIXBASE.TimeInt.fs;
V_care=0.64833*30;
[A_1,A_2,M_new,K_new,C_new,vect_preload]=Get_int_Matrix(FLUTEER_MATRIXBASE,method,V_care);
res=0.72;
B=(1+res)*[0 0 0 C_1 C_2 -1 0 0]';
W=zeros(8);
W(6,6)=B(6);
W(4,6)=B(4);
W(5,6)=B(5);
%%
Scale = 1000000;
%

% InitCond=Static_equilibrium_output(gap,A_1,A_2,vect_preload)
% %

% the start point on the poincare section
x0 = [-0.00465412709916764,0.0606250812329585,-0.00962039197420514,-0.00102284858298862,-0.000663458443753684,0.00466523614885600,0.0624456248892721,-5.96091134697107e-06]';
focus=Static_equilibrium_output(gap,A_1,A_2,vect_preload);
%define the period T (found previously by numerical)
T=0.132671235175032;
%% new definition of state
% Balance_state1=-A_1\vect_preload(:,1);
% correction=Balance_state1;
% correction(3)=-0.01;
% vect_preload(:,1)=A_1*(correction-Balance_state1);
% x0=[3.24890790952170e-06,3.86371182043812e-06,-4.77636050311232e-15,0.00104900482198899,-0.000614240949814067,0.00758675099598211,0.000166291796263342,-3.36393278467710e-05]';
%%
% define the num of loop and store the LPE
delta=0*rand(8,1)/Scale;
InitCond=x0+delta;
N_loop=1;
% for kk=1:N_loop
% tspan=[0 T];
% [tout,yout,yeout0,teout,yeout,ieout,amplitude]=backlash_airfoil_restitution(res,A_1,MHH,InitCond,vect_preload,tspan,gap,fs);
% InitCond=yout(end,:)';
% lpe(kk)=((InitCond-x0)'*delta)/(norm(delta)^2);
% delta=yout(end,:)'-x0;
% end
tspan=[0 N_loop*T];
[tout,yout,yeout0,teout,yeout,ieout,amplitude]=backlash_airfoil_restitution(res,A_1,MHH,InitCond,vect_preload,tspan,gap,fs);
%
X=yout(end,:);
figure
plot(tout,yout(:,3),'r-')
hold on
 if ~isempty(yeout)
plot(teout,yeout(:,3),'ko');
 end
 %
 figure(3)
plot(yout(:,3),yout(:,6))
hold on
plot(InitCond(3),InitCond(6),'r*')
plot(X(3),X(6),'bo')
% Chose a positive definite matrix Q to find corresponding P matrix to
% construct the Lyapunov function 
[V,D]=eig(A_1);
Q=diag(abs(real(diag(D))));
Q=diag([1 2 3 4 5 6 7 8]);
P=lyap(A_1',Q);
% P=inv(V')*inv(V)
Check = A_1*P+P*A_1'
%
x_0=[-0.00463039246834882,0.0606390533238653,-0.0100000000000000,0.000995237839456032,0.000788341253839342,-0.0105389879437051,0.0624496225276524,-3.36552231231508e-05]';
justice=(W*x_0)'*((P+P')*x_0)
% new definition
Composed = W'*(P+P')
diag(eig(Composed))
x_0'*Composed*x_0
focus'*Composed*focus
%% inspect the constructed Lyapunov function along the trajectory
focus=Static_equilibrium_output(gap,A_1,A_2,vect_preload);
Yout=yout;
Lyp=zeros(size(yout,1),1);
Deriv=zeros(size(yout,1),1);
K_Energy=zeros(size(yout,1),1);
P_Energy=zeros(size(yout,1),1);
for kk= 1:size(yout,1)
      Yout(kk,:)=Yout(kk,:)-focus';
    Lyp(kk)=Yout(kk,:)*P*Yout(kk,:)';
    Deriv(kk)=-Yout(kk,:)*Q*Yout(kk,:)';
    K_Energy(kk)= 0.5*yout(kk,4:6)*M_new*yout(kk,4:6)'+0.*yout(kk,1:3)*K_new*yout(kk,1:3)';
    P_Energy(kk)= 0.5*(yout(kk,1:3)-focus(1:3)')*K_new*(yout(kk,1:3)-focus(1:3)')'-0*0.5*focus(1:3)'*K_new*focus(1:3);
end
figure
plot(tout,Lyp)
figure
plot(tout,Deriv)
figure
plot(tout,K_Energy)
figure
plot(tout,P_Energy)
figure
plot(tout,K_Energy+P_Energy)
%
function   equilibrium_state=Static_equilibrium_output(gap,A_1,A_2,vect_preload)
Balance_state1=-A_1\vect_preload(:,1);
    Balance_state2=-A_2\vect_preload(:,1);
    %
    if abs(Balance_state1(2))>gap
        equilibrium_state=-A_2\(vect_preload(:,1)+sign(Balance_state1(2))*vect_preload(:,2));
%             equilibrium_state(2)=gap*sign(Balance_state1(2));
        
    else
        equilibrium_state=Balance_state1;
    end
    %
    equilibrium_state=Balance_state1;
end