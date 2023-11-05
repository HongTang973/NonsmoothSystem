% check if the conjecture is true: p.36 in the book Simpson
%  ---Bifurcation in piecewise-smooth continuous systems ---
% Conjecture :eigenvalue of convex conbination of the Jacobians should
% corssthe imaginary axis when there is bifurcation

% ------  LEFT UNPROVEN -----%
clc
close all
clear
%% CHECK THIS BY AIRFOIL MODEL 
public.plot.gca.linewidth = 3;
public.plot.gca.fontname = 'Times New Roman';
public.plot.gca.fontsize = 16;

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
% V_care=3*30;
% V_care=0.64833*30;
% p1 --point A in ZoomIn 
% V_care=0.6488*30;
% p2 --point B in ZoomIn 
V_care=0.6495*30;

% p3 --point C in ZoomIn 
% V_care=0.6505*30;

% V_care=2*30;
[A_1,A_2,M_new,K_new,C_new,vect_preload]=Get_int_Matrix(FLUTEER_MATRIXBASE,method,V_care);

% FIRST GET THE MATRIX AND THE EIGENVALUES RESPECTIVELY

E1 = eig(A_1);
E2 = eig(A_2);

% 
convex_A = @(s) s*A_2 + (1-s)*A_1;
s_list = 0:0.001:1;
EIG = zeros(length(s_list),size(A_1,1));
EIG(1,:) = E1.';
EIG(end,:) = E2.';
for i = 2:length(s_list)-1
    EIG(i,:)= eig(convex_A(s_list(i))).';
end

% plot 

FIG1 =figure;
plot(real(EIG),imag(EIG),'.','color','blue')
hold on
h1 = plot(real(E1),imag(E1),'o','markersize',4,'markerfacecolor','red','displayname','L');
h2 =plot(real(E2),imag(E2),'o','markersize',4,'markerfacecolor','black','displayname','R');
set(gca,'fontsize',public.plot.gca.fontsize,...
        'fontname',public.plot.gca.fontname,...
        'linewidth',public.plot.gca.linewidth)
legend([h1 h2], 'location','best')
xlabel('Re')
ylabel('Img')
grid on

