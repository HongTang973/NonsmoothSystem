% check if the conjecture is true: p.36 in the book Simpson
%  ---Bifurcation in piecewise-smooth continuous systems ---
% Conjecture :eigenvalue of convex conbination of the Jacobians should
% corssthe imaginary axis when there is bifurcation

% check |exp(A(t)) * exp(B(t)) - I| relationship with |exp(A+B) - I|
% ------  LEFT UNPROVEN -----%
clc
close all
clear
%% CHECK THIS BY AIRFOIL MODEL 
public.plot.gca.linewidth = 3;
public.plot.gca.fontname = 'Times New Roman';
public.plot.gca.fontsize = 16;

%%  define the linear part & get the matrix A
% matrices from CARMONA, et al. 2004

A_1 = [-3.2 -1 0; 25.61 0 -1; -75.03 0 0];
A_2 = [-1 -1 0;1.28 0 -1; -0.624 0 0];


% FIRST GET THE MATRIX AND THE EIGENVALUES RESPECTIVELY

E1 = eig(A_1);
E2 = eig(A_2);

% 
n_m = size(A_1,1);

convex_A = @(s) s*A_2 + (1-s)*A_1;

s_list = 0:0.001:1;

EIG = zeros(length(s_list),n_m);
det_EC = zeros(length(s_list),2);

EIG(1,:) = E1.';
EIG(end,:) = E2.';

for i = 2:length(s_list)-1
    
    temp = convex_A(s_list(i));
    
    EIG(i,:)= eig(temp).';
    
    det_EC(i,1) = det(expm(temp)-eye(n_m));
    
    det_EC(i,2) = det(expm((1-s_list(i))*A_1)*expm(s_list(i)*A_2)-eye(n_m));
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

xlabel('Re')
ylabel('Img')
grid on

% s_list(abs(real(EIG(:,1)))<1e-3)
% s_list(abs(real(EIG(:,2)))<1e-3)
% s_list(abs(real(EIG(:,3)))<1e-3)

% mark the special point with different color
n_ = 3;

len = 1/(n_ + 1);
s_s_list = 0:len:1;

% color 
c_ = (1-s_s_list')*[0 1 0] + s_s_list'*[1 1 0];

% mark the special point
for i = 2:length(s_s_list)-1
    s_EIG= eig(convex_A(s_s_list(i))).';
    plot(real(s_EIG),imag(s_EIG),'o','markersize',4,'markerfacecolor',c_(i,:));
end
s_EIG= eig(convex_A(0.5)).';
plot(real(s_EIG),imag(s_EIG),'o','markersize',6,'markerfacecolor',[0 0 0],'linewidth',2);

legend([h1 h2], 'location','best')


%%
FIG2 = figure;
h3 = plot(s_list,det_EC(:,1),'r-','displayname','Direct Sum');
hold on
h4 = plot(s_list,det_EC(:,2),'k-','displayname','Multiply Sum');
set(gca,'fontsize',public.plot.gca.fontsize,...
        'fontname',public.plot.gca.fontname,...
        'linewidth',public.plot.gca.linewidth)

xlabel('S')
ylabel('det')
grid on
legend([h3 h4], 'location','best')