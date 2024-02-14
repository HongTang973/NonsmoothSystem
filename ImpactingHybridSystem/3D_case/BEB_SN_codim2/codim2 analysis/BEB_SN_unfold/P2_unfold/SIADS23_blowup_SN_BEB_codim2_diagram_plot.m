close all
clc
%
load('./codim2_SN_mu_eta_curve.mat')
filepath = [pwd,'\figures'];
if ~exist(filepath, 'dir')
    mkdir(filepath)
end


BEB_U = 0;
SN_U  = CR_points(2,7);
SN_eta  = CR_points(2,8);
L_U   = -0.04;
R_U   = 0.04;

FIG1 = figure;
%>
Y = [admiss_SN_eta',0.22*ones(size(admiss_SN_mu))'-admiss_SN_eta'];
a = area(admiss_SN_mu,Y,-0.4);
a(1).FaceColor = [1 1 1];
a(1).FaceAlpha = 0.0;
a(2).FaceColor = [0 1 0];
a(2).FaceAlpha = 0.1;
hold on
%
h0 = plot([BEB_U BEB_U],[-0.3 0.2],'b-','LineWidth',2.5,'DisplayName','BEB line');
set(gca,'Layer','top')
h1 = plot(virtual_SN_mu, virtual_SN_eta,'k--','LineWidth',2.5,'DisplayName','Virtual LCO');
h2 = plot(admiss_SN_mu, admiss_SN_eta,'r-','LineWidth',2.5,'DisplayName','Admissible LCO');
box2_x = [L_U L_U BEB_U BEB_U L_U];
box2_y = [-0.3 0.2 0.2 -0.3 -0.3];
V2 = [box2_x(1), box2_y(1);
    box2_x(2), box2_y(2);
    box2_x(3), box2_y(3);
    box2_x(4), box2_y(4)];
F2 = [1 2 3 4];
patch( 'Faces',F2, 'Vertices',V2,'EdgeColor',[1 1 1],'LineStyle','none', 'FaceColor', [1 0 0], 'FaceAlpha', 0.1)
xlim([L_U R_U])
ylim([-0.3 0.2])

grid on
box on



h3 = plot(BEB_U,0,'ko','markersize',6,'markerfacecolor',[0.5 0.5 0.5],'linewidth',2,'DisplayName','BEB point');
plot([-0.01 0.03], [SN_eta SN_eta],'k-','linewidth',2);
h4 = plot(SN_U, SN_eta,'rd','markersize',7,'markerfacecolor',[1 1 1],'linewidth',1.5,'DisplayName','SN point');
plot(0, SN_eta,'ko','markersize',6,'markerfacecolor',[0.5 0.5 0.5],'linewidth',1.5)

xlabel('$\mu$', 'interpreter','latex')
ylabel('$\eta$', 'interpreter','latex')
legend([h0 h1 h2 h3 h4], 'Interpreter','latex','location','best')
adj_plot_theme_I(FIG1)

ax =gca;
ax.gridlinewidth = 1;
% ax.GridLineWidth = 1;
ytick = ax.YTick;
new_label = cell(1,length(ytick));
ytick(4)  = 0;
for i =2:2:length(ytick)
    new_label{i}= sprintf('%g',ytick(i));
end
set(gca,'YTickLabel',new_label);

    name ='3DOF_SIADS23_BEB_SN_codim2';
    picname2 = strcat(filepath,['\',name,'.pdf']);
%    exportgraphics(FIG1,picname2,'ContentType','vector')