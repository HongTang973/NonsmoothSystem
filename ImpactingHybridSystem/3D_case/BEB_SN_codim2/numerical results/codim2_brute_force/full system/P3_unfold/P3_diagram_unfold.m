load('.\P3_brute_force_diagram.mat')

close all

% filepath = [pwd,'\figures'];
% if ~exist(filepath, 'dir')
%     mkdir(filepath)
% end
BEB_U = 0;
SN_U  = 0.0264121728473816;
% SN_eta  = CR_points(2,8);
L_U   = -0.02;
R_U   = 0.03;
L     = -40.404;

%>
FIG1 = figure;
hold on
data             = A_diagram(:,2:end);
data(data<=1e-5) = NaN;
on_eq_id         = A_diagram(:,1)>-SN_U;
% Peter_CreatePlotInOrigin(beta_diagram,'beta_diagram')   
plot(A_diagram(on_eq_id,1)+SN_U,data(on_eq_id),'k.','markersize',12)
h0 = plot(0,0,'k.','linewidth',3,'markersize',10,'displayname','Stable LCO');
hold on

% Color region  to zoom in
box1_x = [BEB_U BEB_U SN_U SN_U BEB_U];
box1_y = [-0.1 0.16 0.16 -0.1 -0.1 ];
V1 = [box1_x(1), box1_y(1);
    box1_x(2), box1_y(2);
    box1_x(3), box1_y(3);
    box1_x(4), box1_y(4)];
F1 = [1 2 3 4];
patch( 'Faces',F1, 'Vertices',V1,'EdgeColor',[1 1 1],'LineStyle','none', 'FaceColor', 'green', 'FaceAlpha', 0.1)

box2_x = [L_U L_U BEB_U BEB_U L_U];
box2_y = [-0.1 0.16 0.16 -0.1 -0.1 ];
V2 = [box2_x(1), box2_y(1);
    box2_x(2), box2_y(2);
    box2_x(3), box2_y(3);
    box2_x(4), box2_y(4)];
F2 = [1 2 3 4];
patch( 'Faces',F2, 'Vertices',V2,'EdgeColor',[1 1 1],'LineStyle','none', 'FaceColor', [1 0 0], 'FaceAlpha', 0.1)
load P3_semi_ana_diagram.mat
S_br = S_amp.*S_u(7,:).*S_B - S_u(7,:).*S_amp ;
h10 = plot(S_u(7,:),S_br,'g-','linewidth', 1., 'displayname','Stable LCO');
hold on
h20 = plot(US_u(7,:),US_amp.*US_u(7,:).*US_B - US_u(7,:).*US_amp,'k--','linewidth', 2,'displayname','Unstable LCO');

h1=plot([0:-0.001:-0.03],L*[0:-0.001:-0.03],'r-','linewidth',3,'Markersize',10,'displayname','Stable AE');
h2=plot([BEB_U R_U],[0 0],'-','color',[0 0 1],'linewidth',3,'Markersize',10,'displayname','Stable PE');
h3 = plot(BEB_U,0,'ko','markersize',6,'markerfacecolor',[0.5 0.5 0.5],'linewidth',2,'DisplayName','BEB point');
h4 = plot(SN_U,S_br(end),'rd','markersize',4,'markerfacecolor',[1 1 1],'linewidth',2,'DisplayName','SN point');
legend([h0 h20 h1 h2 h3 h4], 'Interpreter','latex','location','best')
legend boxoff 

box on
grid on
ax = gca;
ax.GridLineWidth = 1;
xlabel('$\mu$','interpreter','latex');
ylabel('$x_1$','interpreter','latex');
xlim([L_U R_U])
ylim([ -0.01 0.16])
adj_plot_theme_I(FIG1)

xtick = ax.XTick;
x_new_label = cell(1,length(xtick));
for i =2:2:length(xtick)
    x_new_label{i}= sprintf('%g',xtick(i));
end
set(gca,'XTickLabel',x_new_label);
%% post processing
name ='.\3D_BEB_SN_codim2_P2_unfold_full';
    % picname2 = strcat(filepath,['\',name,'.pdf']);
      picname2 = [name,'.pdf'];
%    exportgraphics(FIG1,picname2,'ContentType','vector')