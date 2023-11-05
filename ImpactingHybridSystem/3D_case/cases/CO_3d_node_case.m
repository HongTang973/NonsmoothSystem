% continuation analysis

% at the current parameter point:
% first continuation in one codimension get several candidates of other
% starting parameter points and then do continuation in seperate point


% parameter space : lambda_1, lambda_2, b3

clc
close all
clear 
public.plot.gca.linewidth = 3;
public.plot.gca.fontname = 'Times New Roman';
public.plot.gca.fontsize = 16;
public.plot.gca.ws_x = 0.48;
public.plot.gca.ws_y = 0.3354;
public.plot.gca.width = 4.5208;
public.plot.gca.height = 3.5656;
% define the zero functions
% par = {'lambda_1', 'lambda_2', 'b2','b3', 'T'}
% 
%  par =  [-1, 0.2, -1, 2.5, 1, 0.875926384431162];
 par = [-11.8818;0.2;-1;2.5;1;0.1];
% par1 =  [-1, 0.2, -1,2.5, 1, 0.7627];
% par2 =  [-1, 0.001, 2.5, 1, 0.7627];
ind_p =1; % the lambda_1
ind_x =2; % the Period T
dir =1;
% [prob_det] = ZeroFunctions_3D_node([-1, -0.01, 2.5, 1, 0.7627])
p_span = [-12  12];

x_span = [-12 12];

tic
par_1 = [-11.8818;0.2;-1;2.5;1;0.1];
p_span_1 = [-20  -0.01];
x_span_1 = [0.1 10];
[u1,iter1] = codim1_PC(@ZeroFunctions_3D_IHS,par_1,ind_p,ind_x,p_span_1,x_span_1,-1);
FIG1 = figure;set(FIG1,'Units','inches');
hax1 = axes;
plot(hax1,u1(ind_p,:),u1(ind_x,:),'.','color','black','linewidth',1.5)
hold on
plot(hax1,u1(ind_p,1),u1(ind_x,1),'o','color','red','linewidth',2,'markersize',8,'markerfacecolor',[1 1 1])
%
par_2 = [-11.8818;0.2;-1;2.5;1;0.1];
p_span_2 = [-12  -0.01];
x_span_2 = [0.01 12];
[u2,iter2] = codim1_PC(@ZeroFunctions_3D_IHS,par_2,ind_p,ind_x,p_span_2,x_span_2,1);
plot(hax1,u2(ind_p,:),u2(ind_x,:),'.','color','green','linewidth',1.5)
%
par_3 = [-11.6818;-0.00727;-1;2.5;1;0.1];
p_span_3 = [-12  -0.01];
x_span_3 = [-12 -0.001];
[u3,iter3] = codim1_PC(@ZeroFunctions_3D_IHS,par_3,ind_p,ind_x,p_span_3,x_span_3,-1);
plot(hax1,u3(ind_p,:),u3(ind_x,:),'.','color','blue','linewidth',1.5)
%
par_4 = [0.009105;-11.66;-1;2.5;1;0.1];
p_span_4 = [0.0  12];
x_span_4 = [-22 -0.001];
[u4,iter4] = codim1_PC(@ZeroFunctions_3D_IHS,par_4,ind_p,ind_x,p_span_4,x_span_4,1);
plot(hax1,u4(ind_p,:),u4(ind_x,:),'.','color','black','linewidth',1.5)

xlabel('$\lambda_2$','interpreter','latex')
ylabel('$T$','interpreter','latex')

set(gca,'fontsize',public.plot.gca.fontsize,...
        'fontname',public.plot.gca.fontname,...
        'linewidth',public.plot.gca.linewidth,...
        'Units','inches')
% [u1,iter1] = codim1_PC_skip0(@ZeroFunctions_3D_IHS,par,ind_p,ind_x,p_span,x_span,-1);
% [u2,iter2] = codim1_PC_skip0(@ZeroFunctions_3D_IHS,par,ind_p,ind_x,p_span,x_span,1);
% [u3,iter3] = codim1_PC(@ZeroFunctions_3D_IHS,par1,ind_p,ind_x,p_span,x_span,-1);

CPU_tim = toc

u_all = [u1,u2,u3,u4];
equi_type =1;
bifur_index = [];
FloquetM = [];
Sign_V =[];
for i = 1:size(u_all,2)
    par = u_all(:,i);
    
    lambda_1 = par(1);
    
    lambda_2 = par(2);
    
    lambda_3 = par(3);
    
    a1 = lambda_1 + lambda_2 + lambda_3;
    a2 = -(lambda_1 * lambda_2 + lambda_2*lambda_3 + lambda_1*lambda_3);
    a3 = lambda_1 * lambda_2 * lambda_3;
    %
    b2 = par(4);
    b3 =par(5);
    
    A = [a1 1 0;a2 0 1; a3 0 0];
    
    
    B = [0; b2; b3];
    
    
    C = [1,0,0];
    
    R = eye(length(C)) - B*C*A;
    
    T = par(6);
    
   [~,~,~,LCO,~] = LCO_Det_search(T,R,A,C,equi_type);
    [Mono_p,Salt_p]=Floque_Multipliers(T,LCO,R,A,C) ;  
    
        if max(abs(Salt_p)) >1+0.001%  && max(abs(Salt_p)) -1 < 1e-6
           bifur_index = [bifur_index,i];
            FloquetM = [FloquetM,sort(Salt_p,'descend')];
        end
    Sign_V = [Sign_V,sign(C*A*LCO)];
end
u_all(:,2461)
u_all(:,2460)
%% step 2: expand on another codimension: seven branches
ind_p =2; % the lambda_2
ind_x =5; % the b_3

%% 
par_1 = branch_1_0;

p_span_1_1 = [0.01 10];
x_span_1_1 = [-2 5];
[u_2_b3_1_1,iter1] = codim1_PC(@ZeroFunctions_3D_IHS,par_1,ind_p,ind_x,p_span_1_1,x_span_1_1,-1);
FIG2 = figure;set(FIG2,'Units','inches');
hax2 = axes;
plot(hax2,u_2_b3_1_1(ind_p,:),u_2_b3_1_1(ind_x,:),'.','color','black','linewidth',1.5)
hold on
plot(hax2,u_2_b3_1_1(ind_p,1),u_2_b3_1_1(ind_x,1),'o','color','blue','linewidth',2,'markersize',8,'markerfacecolor',[1 1 1])

p_span_1_2 = [0.01 10];
x_span_1_2 = [-2 5];

[u_2_b3_1_2,iter1] = codim1_PC(@ZeroFunctions_3D_IHS,par_1,ind_p,ind_x,p_span_1_1,x_span_1_1,1);

plot(hax2,u_2_b3_1_2(ind_p,:),u_2_b3_1_2(ind_x,:),'.','color','black','linewidth',1.5)

u_2_b3_1 = [u_2_b3_1_1,u_2_b3_1_2];

%%  *****
par_2 = branch_2_0;

p_span_2_1 = [0.01 10];
x_span_2_1 = [-2 5];
[u_2_b3_2_1,iter1] = codim1_PC(@ZeroFunctions_3D_IHS,par_2,ind_p,ind_x,p_span_2_1,x_span_2_1,-1);

plot(hax2,u_2_b3_2_1(ind_p,:),u_2_b3_2_1(ind_x,:),'.','color','black','linewidth',1.5)
hold on
plot(hax2,u_2_b3_2_1(ind_p,1),u_2_b3_2_1(ind_x,1),'o','color','blue','linewidth',2,'markersize',8,'markerfacecolor',[1 1 1])

p_span_2_2 = [0.01 10];
x_span_2_2 = [-2 5];

[u_2_b3_2_2,iter1] = codim1_PC(@ZeroFunctions_3D_IHS,par_2,ind_p,ind_x,p_span_2_1,x_span_2_1,1);

plot(hax2,u_2_b3_2_2(ind_p,:),u_2_b3_2_2(ind_x,:),'.','color','black','linewidth',1.5)

u_2_b3_2 = [u_2_b3_2_1,u_2_b3_2_2];


%

par_3 = branch_3_0;

p_span_3_1 = [0.01 10];
x_span_3_1 = [-2 5];
[u_2_b3_3_1,iter1] = codim1_PC(@ZeroFunctions_3D_IHS,par_3,ind_p,ind_x,p_span_3_1,x_span_3_1,-1);

plot(hax2,u_2_b3_3_1(ind_p,:),u_2_b3_3_1(ind_x,:),'.','color','black','linewidth',1.5)
hold on
plot(hax2,u_2_b3_3_1(ind_p,1),u_2_b3_3_1(ind_x,1),'o','color','blue','linewidth',2,'markersize',8,'markerfacecolor',[1 1 1])

p_span_3_2 = [0.01 10];
x_span_3_2 = [-2 5];

[u_2_b3_3_2,iter1] = codim1_PC(@ZeroFunctions_3D_IHS,par_3,ind_p,ind_x,p_span_3_1,x_span_3_1,1);

plot(hax2,u_2_b3_3_2(ind_p,:),u_2_b3_3_2(ind_x,:),'.','color','black','linewidth',1.5)

%
u_2_b3_3 = [u_2_b3_3_1,u_2_b3_3_2];
%
par_4 = branch_4_0;

p_span_4_1 = [0.01 10];
x_span_4_1 = [-2 5];
[u_2_b3_4_1,iter1] = codim1_PC(@ZeroFunctions_3D_IHS,par_4,ind_p,ind_x,p_span_4_1,x_span_4_1,-1);

plot(hax2,u_2_b3_4_1(ind_p,:),u_2_b3_4_1(ind_x,:),'.','color','black','linewidth',1.5)
hold on
plot(hax2,u_2_b3_4_1(ind_p,1),u_2_b3_4_1(ind_x,1),'o','color','blue','linewidth',2,'markersize',8,'markerfacecolor',[1 1 1])

p_span_3_2 = [0.01 10];
x_span_3_2 = [-2 5];

[u_2_b3_4_2,iter1] = codim1_PC(@ZeroFunctions_3D_IHS,par_4,ind_p,ind_x,p_span_4_1,x_span_4_1,1);

plot(hax2,u_2_b3_4_2(ind_p,:),u_2_b3_4_2(ind_x,:),'.','color','black','linewidth',1.5)

u_2_b3_4 = [u_2_b3_4_1,u_2_b3_4_2];

%
par_5 = branch_5_0;

p_span_5_1 = [-4 -0.01];
x_span_5_1 = [-2 5];
[u_2_b3_5_1,iter1] = codim1_PC(@ZeroFunctions_3D_IHS,par_5,ind_p,ind_x,p_span_5_1,x_span_5_1,-1);

plot(hax2,u_2_b3_5_1(ind_p,:),u_2_b3_5_1(ind_x,:),'.','color','black','linewidth',1.5)
hold on
plot(hax2,u_2_b3_5_1(ind_p,1),u_2_b3_5_1(ind_x,1),'o','color','blue','linewidth',2,'markersize',8,'markerfacecolor',[1 1 1])

p_span_5_2 = [0.001 10];
x_span_5_2 = [-2 5];

[u_2_b3_5_2,iter1] = codim1_PC(@ZeroFunctions_3D_IHS,[-11.6456;0.01;-1.0000;2.5000;0.9402;0.1000],ind_p,ind_x,p_span_5_2,x_span_5_2,1);

plot(hax2,u_2_b3_5_2(ind_p,:),u_2_b3_5_2(ind_x,:),'.','color','black','linewidth',1.5)


p_span_5_3 = [-4 -0.01];
x_span_5_3 = [-2 5];
[u_2_b3_5_3,iter1] = codim1_PC(@ZeroFunctions_3D_IHS,par_5,ind_p,ind_x,p_span_5_1,x_span_5_1,1);

plot(hax2,u_2_b3_5_3(ind_p,:),u_2_b3_5_3(ind_x,:),'.','color','black','linewidth',1.5)

u_2_b3_5 = [u_2_b3_5_1,u_2_b3_5_2,u_2_b3_5_2];
%

%
par_6 = branch_6_0;

p_span_6_1 = [-4 -0.01];
x_span_6_1 = [-2 5];
[u_2_b3_6_1,iter1] = codim1_PC(@ZeroFunctions_3D_IHS,par_6,ind_p,ind_x,p_span_6_1,x_span_6_1,-1);

plot(hax2,u_2_b3_6_1(ind_p,:),u_2_b3_6_1(ind_x,:),'.','color','black','linewidth',1.5)
hold on
plot(hax2,u_2_b3_6_1(ind_p,1),u_2_b3_6_1(ind_x,1),'o','color','blue','linewidth',2,'markersize',8,'markerfacecolor',[1 1 1])
%
p_span_6_2 = [-4 10];
x_span_6_2 = [-2 5];

[u_2_b3_6_2,iter1] = codim1_PC(@ZeroFunctions_3D_IHS,par_6,ind_p,ind_x,p_span_6_2,x_span_6_2,1);

plot(hax2,u_2_b3_6_2(ind_p,:),u_2_b3_6_2(ind_x,:),'.','color','black','linewidth',1.5)
%
p_span_6_3 = [-4 10];
x_span_6_3 = [-2 5];

[u_2_b3_6_3,iter1] = codim1_PC(@ZeroFunctions_3D_IHS,[-11.5078;0.004637;-1.0000;2.5000;0.5984;0.1000],ind_p,ind_x,p_span_6_3,x_span_6_3,1);

plot(hax2,u_2_b3_6_3(ind_p,:),u_2_b3_6_3(ind_x,:),'.','color','black','linewidth',1.5)

u_2_b3_6 = [u_2_b3_6_1,u_2_b3_6_2,u_2_b3_6_3];
%
%
par_7 = branch_7_0;

p_span_7_1 = [-4 -0.01];
x_span_7_1 = [-2 5];
[u_2_b3_7_1,iter1] = codim1_PC(@ZeroFunctions_3D_IHS,par_7,ind_p,ind_x,p_span_7_1,x_span_7_1,-1);

plot(hax2,u_2_b3_7_1(ind_p,:),u_2_b3_7_1(ind_x,:),'.','color','black','linewidth',1.5)
hold on
plot(hax2,u_2_b3_7_1(ind_p,1),u_2_b3_7_1(ind_x,1),'o','color','blue','linewidth',2,'markersize',8,'markerfacecolor',[1 1 1])

%

p_span_7_2 = [-4 -0.01];
x_span_7_2 = [-2 5];

[u_2_b3_7_2,iter1] = codim1_PC(@ZeroFunctions_3D_IHS,par_7,ind_p,ind_x,p_span_7_2,x_span_7_2,1);

plot(hax2,u_2_b3_7_2(ind_p,:),u_2_b3_7_2(ind_x,:),'.','color','black','linewidth',1.5)

%
p_span_7_3 = [0.01 10];
x_span_7_3 = [-2 5];

[u_2_b3_7_3,iter1] = codim1_PC(@ZeroFunctions_3D_IHS,[-11.3560144011424;0.02;-1;2.5;0.15;0.1],ind_p,ind_x,p_span_7_3,x_span_7_3,1);

plot(hax2,u_2_b3_7_3(ind_p,:),u_2_b3_7_3(ind_x,:),'.','color','black','linewidth',1.5)

u_2_b3_7 = [u_2_b3_7_1,u_2_b3_7_2,u_2_b3_7_3];
%


figure
plot3(u_all(1,:),u_all(2,:),u_all(5,:),'.')
hold on
plot3(u_2_b3_1(1,:),u_2_b3_1(2,:),u_2_b3_1(5,:),'.')
plot3(u_2_b3_2(1,:),u_2_b3_2(2,:),u_2_b3_2(5,:),'.')
plot3(u_2_b3_3(1,:),u_2_b3_3(2,:),u_2_b3_3(5,:),'.')
plot3(u_2_b3_4(1,:),u_2_b3_4(2,:),u_2_b3_4(5,:),'.')
plot3(u_2_b3_5(1,:),u_2_b3_5(2,:),u_2_b3_5(5,:),'.')
plot3(u_2_b3_6(1,:),u_2_b3_6(2,:),u_2_b3_6(5,:),'.')
plot3(u_2_b3_7(1,:),u_2_b3_7(2,:),u_2_b3_7(5,:),'.')
xlabel('$\lambda_1$','interpreter','latex')
ylabel('$\lambda_2$','interpreter','latex')
zlabel('$b_3$','interpreter','latex')
xlim([-14 -10])
ylim([-2.5 2.5])
zlim([-1 2.5])


