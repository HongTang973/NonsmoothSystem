% continuation analysis

% at the current parameter point:
% first continuation in one codimension get several candidates of other
% starting parameter points and then do continuation in seperate point


% parameter space : lambda_1, lambda_2, b3

clc
close all

% define the zero functions

% 
% par = {'lambda_1', 'lambda_2','lambda_3', 'b2','b3', 'T'}
% par =  [-1, 0.2, -1, 2.5, 0.5, 0.875926384431162]; 
% initial point from algorithm 

%% First step: vary lambda_1 tp get the desired T = epsilon_2

% continuation to a new point with: T = epsilon_2 = 0.1;

par = [-11.8818;0.2;-1;2.5;1;0.1];



%% Second step: fix the period, vary the lambda_1 to get some 
% lambda_1/lambda_2 curve 
ind_p =2; % the lambda_2
ind_x =5; % the Period T
dir =1;
% [prob_det] = ZeroFunctions_3D_node([-1, -0.01, 2.5, 1, 0.7627])
p_span = [-10  1];

x_span = [0 10];

tic
[u1,iter1] = codim1_PC(@ZeroFunctions_3D_IHS,par,ind_p,ind_x,p_span,x_span,-1);
[u2,iter2] = codim1_PC(@ZeroFunctions_3D_IHS,par,ind_p,ind_x,p_span,x_span,1);
[u3,iter3] = codim1_PC(@ZeroFunctions_3D_IHS,par,ind_p,ind_x,p_span,x_span,-1);

CPU_tim = toc

FIG1 = figure;set(FIG1,'Units','inches');
hax1 = axes;
plot(u1(ind_p,:),u1(end,:),'color','black','linewidth',1.5)
hold on
plot(u2(ind_p,:),u2(end,:),'color','blue','linewidth',1.5)
plot(u3(ind_p,:),u3(end,:),'color','blue','linewidth',1.5)
plot(u1(ind_p,1),u1(end,1),'o','color','red','linewidth',2,'markersize',8,'markerfacecolor',[1 1 1])
xlabel('$\lambda_2$','interpreter','latex')
ylabel('$T$','interpreter','latex')

set(gca,'fontsize',public.plot.gca.fontsize,...
        'fontname',public.plot.gca.fontname,...
        'linewidth',public.plot.gca.linewidth,...
        'Units','inches')