clc
% clear
close all
load('.\BEB_2_SN_LCO_continuation_P3.mat')
% P3 = [-0.1	0.2	-0.5	1.78192697901049	1.6	4.7959612404933 0.0264121728473816  0.1];
P3 = [-0.1,0.2,-0.5,1.78192697901049,1.6,4.73578626452264,0.0264121728473816,0.1];
mu_crit      = P3(7);
eta_crit     = P3(8);
SN_point_admis = P3;
% mu_crit  = 0.0264121728473816;
% eta_crit = 0.1;
SN_point_admis(3:4) = SN_point_admis(3:4) + [-mu_crit,eta_crit];
% mu_crit = 0.0264121728473816;
% SN_point_admis = [-0.1  0.2	-0.5 1.78192697901049 	1.6	5.84686867166138 ];
%> use the numerics to evaluate the local expansions
% SN_SIADS23_norm_form;


%% > observe the evolution of the fixed point
S_mu_list               = S_u(7,:);
US_mu_list              = US_u(7,:);
S_T_list                = S_u(6,:);
US_T_list               = US_u(6,:);
%>
S_n_num                 = length(S_mu_list);
US_n_num                = length(US_mu_list);
S_par_matrix            = repmat(P3,S_n_num,1);
US_par_matrix           = repmat(P3,US_n_num,1);
S_par_matrix(:,6)       = S_T_list;
S_par_matrix(:,7)       = S_mu_list;
US_par_matrix(:,6)      = US_T_list;
US_par_matrix(:,7)      = US_mu_list;
%> get the fixed point 
[A,B,C,R]   = par2NForm_DummyVar(P3);
% par_matrix(31,:); %> where the irregular point emerges
 x_00       = IC_generator(P3(6),R,A,C,1);
%  u_p        = x_00(2:3);
[Mono_p,Salt_p]     = IC2Floque_Multipliers(P3(6),x_00,R,A,C);
x_0         = R\x_00;
L0          = C*A*x_0;
%>
u_p         = x_0(2:3);

%> added on 16th/Dec/2024
% d      = 0;
% a1 =0;
% c_     = @(mu) sqrt(-3*a1*d*(mu_crit - mu)*0 + c^2);
on_off = 1;
% z_curve_para_1 = @(mu) sqrt(-a0./c_(mu).*(mu_crit - mu)   + ...
%     1*on_off/4./c_(mu).^2 .*(0*a1^2.*(mu_crit - mu).^2 + 0*4*a1.*(mu_crit - mu))) - on_off*a1/2./c_(mu).*(mu_crit - mu);
% z_curve_para_2 = @(mu) -sqrt( -a0./c_(mu).*(mu_crit - mu) + ...
%     1*on_off/4./c_(mu).^2 .*( 0*a1^2.*(mu_crit - mu).^2 + 0*4*a1.*(mu_crit - mu))) - on_off*a1/2./c_(mu).*(mu_crit - mu);
z_curve_para_1 = @(mu) sqrt(a0*c.*(mu - mu_crit) - 0*a1^2*(mu-mu_crit).^2)/c -a1/2./c.*(mu_crit - mu);
z_curve_para_2 = @(mu) -sqrt(a0*c.*(mu - mu_crit)- 0*a1^2*(mu-mu_crit).^2)/c -a1/2./c.*(mu_crit - mu);

US_predict_z = z_curve_para_1(US_mu_list);
S_predict_z  = z_curve_para_2(S_mu_list);
%
S_hat_u  = [];
S_proj_u = [];
S_norm_u = [];
S_v      = [];
S_v_pre  = [];
GM     = [];
for i =1 : S_n_num
    [A,B,C,R]           = par2NForm_DummyVar(S_par_matrix(i,:));
    x_00                = IC_generator(S_par_matrix(i,6),R,A,C,1);
%     [Mono_p,Salt_p]     = IC2Floque_Multipliers(S_par_matrix(i,6),x_00,R,A,C);
%     GM                  = [GM, max(abs(Salt_p))];
%     hat_u       = [hat_u, x_00(2:3)];
%     norm_u      = [norm_u,norm(x_00(2:3) - u_p)];
%     proj_u      = [proj_u,w'*(x_00(2:3) - u_p)];
    x_0         = R\x_00;
    S_hat_u       = [S_hat_u,   x_0(2:3)];
    S_norm_u      = [S_norm_u,  norm(x_0(2:3) - u_p)];
    S_proj_u      = [S_proj_u,  w'*(x_0(2:3) - u_p)];
    S_v           = [S_v,  C*A*x_0];
    S_v_pre           = [S_v_pre,  C*A*[1;u_p+S_predict_z(i)*v]];
end

US_hat_u  = [];
US_proj_u = [];
US_norm_u = [];
US_v      = [];
US_v_pre  = [];
for i =1 : US_n_num
    [A,B,C,R]           = par2NForm_DummyVar(US_par_matrix(i,:));
    x_00                = IC_generator(US_par_matrix(i,6),R,A,C,1);
%     hat_u       = [hat_u, x_00(2:3)];
%     norm_u      = [norm_u,norm(x_00(2:3) - u_p)];
%     proj_u      = [proj_u,w'*(x_00(2:3) - u_p)];
    x_0         = R\x_00;
    US_hat_u       = [US_hat_u,   x_0(2:3)];
    US_norm_u      = [US_norm_u,  norm(x_0(2:3) - u_p)];
    US_proj_u      = [US_proj_u,  w'*(x_0(2:3) - u_p)];
    US_v           = [US_v,  C*A*x_0];
    US_v_pre           = [US_v_pre,  C*A*[1;u_p+US_predict_z(i)*v]];
end

%>

%>  plot the distribution of the fixed points on the Sigma
figure;clf
plot(S_hat_u(1,:), S_hat_u(2,:),'b.')
hold on
plot(US_hat_u(1,:), US_hat_u(2,:),'r.')
plot(u_p(1), u_p(2), 'ro')
plot([u_p(1),u_p(1)+0.02*v(1)],[u_p(2),u_p(2)+0.02*v(2)],'k-' )
plot([u_p(1),u_p(1)-0.02*v(1)],[u_p(2),u_p(2)-0.02*v(2)],'k-' )
title('The fixed points around SN point')
%> 
mu_at_FP                = real(S_mu_list(end));
norm_at_FP              = real(S_norm_u(end));
%>
FIG0 =figure;clf
h1 = plot(S_mu_list,  S_v,'k-','LineWidth',2,'displayname','Stable LCO');
hold on
h2 = plot(US_mu_list, US_v,'k--','LineWidth',2,'displayname','Unstable LCO');
h3 = plot(S_mu_list,  S_v_pre,'r-','LineWidth',2,'displayname','Stable LCO');
h4 = plot(US_mu_list, US_v_pre,'r--','LineWidth',2,'displayname','Unstable LCO');
h5 = plot(mu_at_FP, L0, 'rd','MarkerSize',6,'MarkerFaceColor',[1 1 1],'linewidth',1.5,'displayname','SN');
xlabel('$\mu$','Interpreter','latex')
ylabel('$\hat{v}^-$','Interpreter','latex')
legend([h1 h2 h3 h4 h5],'location','best')
grid on
xlim([0 0.03])
ax = gca;
ax.GridLineWidth = 1;
ytick = ax.YTick;
expVal = -2;
new_label = cell(1,length(ytick));
for i =2:2:length(ytick)
    new_label{i}= sprintf('%g',ytick(i)/(10^expVal));
end
set(gca,'YTickLabel',new_label);
% offset = 0.00;  %how far to the right you want the exponent
annotation('textbox',[0.00 0.85, 0.2, 0.2],...
    'Units', 'pixels',...
    'String',['$\times 10^{' num2str(expVal) '}$'],...
    'fontsize',16,...
    'fontname','times new roman',...
    'Interpreter','latex',...
    'VerticalAlignment','bottom',...
    'EdgeColor','none')

adj_plot_theme_I(FIG0)
%exportgraphics(FIG0,'./Codim2_3D_unfold_SN_BEB_blow_up.pdf','ContentType','vector')

FIG1 = figure;clf
hold on
plot(S_mu_list, S_proj_u, 'b.')
plot(US_mu_list, US_proj_u, 'r.')
plot(mu_at_FP, norm_at_FP, 'ro')
% plot(mu_list_bd0, z_curve_para(mu_list_bd0), 'k--')
% plot(mu_list_bd0, -z_curve_para(mu_list_bd0), 'k--')
plot(US_mu_list, z_curve_para_1(US_mu_list), 'g--')
plot(S_mu_list, z_curve_para_2(S_mu_list), 'k-')
xlabel('$\mu$','Interpreter','latex')
ylabel('$w^{\top}(\hat{u} - u_p)$','Interpreter','latex')
adj_plot_theme_I(FIG1)


%% Approximate the velocity 
% load P3_semi_ana_diagram.mat
L2      = C*inv(A)*M;
SN_mu   = mu_crit;
R_mu    = 0.05;
L_mu    = -0.005;
FIG2 = figure;
h1 = plot(S_mu_list,S_mu_list.*S_v_pre*L2,'r-','LineWidth',2,'displayname','Stable LCO');
hold on

% Color region  to zoom in
box1_x = [SN_mu SN_mu R_mu R_mu SN_mu];
box1_y = [-2.5 0.2 0.2 -2.5 -2.5 ]*1e-3*L2;
V = [box1_x(1), box1_y(1);
    box1_x(2), box1_y(2);
    box1_x(3), box1_y(3);
    box1_x(4), box1_y(4)];
F = [1 2 3 4];
patch( 'Faces',F, 'Vertices',V,'EdgeColor',[1 1 1],'LineStyle','none', 'FaceColor', 'green', 'FaceAlpha', 0.1)

% Color region  to zoom in
box2_x = [L_mu L_mu 0 0 L_mu];
box2_y = [-2.5 0.2 0.2 -2.5 -2.5 ]*1e-3*L2;
V2 = [box2_x(1), box2_y(1);
    box2_x(2), box2_y(2);
    box2_x(3), box2_y(3);
    box2_x(4), box2_y(4)];
F2 = [1 2 3 4];
patch( 'Faces',F2, 'Vertices',V2,'EdgeColor',[1 1 1],'LineStyle','none', 'FaceColor', 'red', 'FaceAlpha', 0.1)

%
h2 = plot(S_mu_list,S_mu_list.*S_v*L2,'k-','LineWidth',2,'displayname','Stable LCO');
h3 = plot(US_mu_list,US_mu_list.*US_v_pre*L2,'r--','LineWidth',2,'displayname','Unstable LCO');
h4 = plot(US_mu_list,US_mu_list.*US_v*L2,'k--','LineWidth',2,'displayname','Unstable LCO');
h5 = plot(mu_crit, mu_crit.*US_v(end)*L2, 'rd','MarkerSize',6,'MarkerFaceColor',[1 1 1],'linewidth',1.5,'displayname','SN');
h6 = plot([0 2*mu_crit], [0 0], 'b-','LineWidth',3,'displayname','Stable PE');
h7 = plot(0,0,'ko','markersize',6,'markerfacecolor',[0.5 0.5 0.5],'linewidth',1.5,'displayname','BEB');
legend([h1 h2 h3 h4 h5 h6 h7],'location','best')
grid on
xlim([-0.005 0.05])
ylim([-2.5 0.1]*1e-3*L2)
xlabel('$\mu$','Interpreter','latex')
ylabel('$v^-$','Interpreter','latex')
adj_plot_theme_I(FIG2)
ax =gca;
ax.GridLineWidth = 0.8;
ytick = ax.YTick;
expVal = -3;
new_label = cell(1,length(ytick));
for i =2:2:length(ytick)
    new_label{i}= sprintf('%g',ytick(i)/(10^expVal));
end
set(gca,'YTickLabel',new_label);
% offset = 0.00;  %how far to the right you want the exponent
annotation('textbox',[0.00 0.85, 0.2, 0.2],...
    'Units', 'pixels',...
    'String',['$\times 10^{' num2str(expVal) '}$'],...
    'fontsize',16,...
    'fontname','times new roman',...
    'Interpreter','latex',...
    'VerticalAlignment','bottom',...
    'EdgeColor','none')
% exportgraphics(FIG2,'./Codim2_3D_unfold_SN_BEB_OP.pdf','ContentType','vector')
%> illustrate the projection 
% FIG2 = figure;
% % L0   = L0 - 1;
% US_mu_shift  = sort(US_mu_list,'descend');
% S_mu_shift  = sort(S_mu_list,'descend');
% L1      = v(1);
% US_v    = US_mu_list.*( (-0.7 - US_mu_shift)+ u_p(1) + L1*US_proj_u);
% S_v     = S_mu_list.*( (-0.7 - S_mu_shift)+ u_p(1) +L1*S_proj_u);
% plot(S_mu_list,S_v)
% hold on
% plot(US_mu_list, US_v)
% 
% h1 = plot(US_mu_list, L0*US_mu_list+US_mu_shift.*US_mu_list + L1*US_mu_list.*US_proj_u,'b--','LineWidth',1.4,'displayname','Unstable LCO');
% hold on
% h2 = plot(S_mu_list, L0*S_mu_list -S_mu_shift.*S_mu_list+ L1*S_mu_list.*S_proj_u,'b-','LineWidth',1.4,'displayname','Stable LCO');
% h3 = plot(US_mu_list, L0*US_mu_list+US_mu_shift.*US_mu_list + L1*US_mu_list.*z_curve_para_1(US_mu_list), 'r--','LineWidth',1.4,'displayname','Unstable LCO');
% h4 = plot(S_mu_list, L0*S_mu_list -S_mu_shift.*S_mu_list+ L1*S_mu_list.*z_curve_para_2(S_mu_list), 'r-','LineWidth',1.4,'displayname','Stable LCO');
% h5 = plot(mu_crit, L0*mu_crit + L1*mu_crit.*z_curve_para_2(mu_crit), 'kd','MarkerSize',5,'MarkerFaceColor',[1 1 1],'displayname','SN');
% h6 = plot(S_mu_list, 0*S_mu_list, 'g-','LineWidth',2,'displayname','Stable PE');
% h7 = plot(-S_mu_list, 0*S_mu_list, 'y-','LineWidth',2,'displayname','Stable AE');
% h8 = plot(0,0,'ko','MarkerSize',5,'MarkerFaceColor',[1 1 1],'displayname','BEB');
% legend([h1 h2 h3 h4 h5 h6 h7,h8],'location','best')
% grid on
% xlim([-0.005 0.03])
% ylim([-2.5 0.1]*1e-3)
% xlabel('$\mu$','Interpreter','latex')
% ylabel('$\mathcal{A}$','Interpreter','latex')
% adj_plot_theme_I(FIG2)
% % exportgraphics(FIG2,'./Codim2_3D_unfold_OP.pdf','ContentType','vector')
FIG3 = figure;
h1 = plot(S_mu_list, S_hat_u(1,:) - u_p(1,:), 'b--','LineWidth',1.4,'displayname','Stable LCO');
hold on
h2 = plot(US_mu_list, US_hat_u(1,:)-u_p(1,:), 'b-','LineWidth',1.4,'displayname','Unstable LCO');

h3 = plot(S_mu_list,  v(1)*z_curve_para_1(S_mu_list), 'r--','LineWidth',1.4,'displayname','Stable LCO');
h4 = plot(US_mu_list, v(1)*z_curve_para_2(US_mu_list), 'r-','LineWidth',1.4,'displayname','Unstable LCO');

%>
h5 = plot(mu_crit, 0, 'ko','MarkerSize',5,'MarkerFaceColor',[1 1 1],'displayname','SN');

xlabel('$\mu$','Interpreter','latex')
ylabel('$\delta(\hat{\mathcal{A}})$','Interpreter','latex')
legend([h1 h2 h3 h4 h5],'location','best')
grid on
ax = gca;
ax.GridLineWidth = 1;
adj_plot_theme_I(FIG3)
% exportgraphics(FIG3,'./Codim2_3D_unfold_scaled.pdf','ContentType','vector')
