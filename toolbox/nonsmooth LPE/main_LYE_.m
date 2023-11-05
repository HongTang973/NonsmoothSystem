% Get the Lyapunov exponent of the Composed map 


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
fs=FLUTEER_MATRIXBASE.TimeInt.fs;
V_care=0.64833*30;
[A_1,A_2,M_new,K_new,C_new,vect_preload]=Get_int_Matrix(FLUTEER_MATRIXBASE,method,V_care);
res=0.72;
%%
Scale = 1000000;
%

% the start point on the poincare section
x0 = [-0.00465412709916764,0.0606250812329585,-0.00962039197420514,-0.00102284858298862,-0.000663458443753684,0.00466523614885600,0.0624456248892721,-5.96091134697107e-06]';
 
%define the period T (found previously by numerical)
T=0.132671235175032;

% define the num of loop and store the LPE
delta=rand(8,1)/Scale;
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
[tout,yout,yeout0,teout,yeout,ieout,amplitude]=backlash_airfoil_restitution(res,A_1,MHH,InitCond,vect_preload,tspan,gap,fs,1);
%
figure
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
X=yout(end,:);
%
figure(3)
plot(yout(:,3),yout(:,6))
hold on
plot(InitCond(3),InitCond(6),'r*')
plot(X(3),X(6),'bo')

%%
(norm(X-x0)/norm(delta))^(1/N_loop)