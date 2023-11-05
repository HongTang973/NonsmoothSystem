% after the line searching code, run this code to get the delicate plot of
% indicating function

close all

rootpath = erase(mfilename('fullpath'),mfilename());
filepath = [rootpath,'\figures'];
if ~exist(filepath, 'dir')
       mkdir(filepath)
end

public.plot.gca.linewidth = 3;
public.plot.gca.fontname = 'Times New Roman';
public.plot.gca.fontsize = 18;
public.plot.gca.ws_x = 0.9;
public.plot.gca.ws_y = 0.7;
public.plot.gca.width = 4.5208;
public.plot.gca.height = 3.5656;
%
FIG1 = figure;set(FIG1,'units','inches')
h1=plot( T,F_1,'r-','linewidth',2,'displayname','$p(t)$');
hold on
% plot( T/pi,V_sign,'-','linewidth',1.2,'displayname','Sign of Velocity')
h2=plot([a b]*pi,[0 0],'--','color',[0 0 0],'linewidth',1.2,'displayname','zero line');
h3=plot(T_chosen(1),F1_chosen(1),'bo','linewidth',2,'markersize',10,...
    'markerfacecolor',[1 1 1],'displayname',['the ',num2str(1),'st candidate']);
h4=plot(T_chosen(2),F1_chosen(2),'bo','markersize',10,'markerfacecolor',...
    'blue','displayname',['the ',num2str(2),'nd candidate']);


xlim([0.,3]*pi)
ylim([-0.02,0.02])
%
xlabel('$t$','interpreter','latex')
ylabel('$p(t)$','interpreter','latex')



legend([h1 h3 h4],'Location', 'best','interpreter','latex')
set(gca,'fontsize',public.plot.gca.fontsize,...
        'fontname',public.plot.gca.fontname,...
        'linewidth',public.plot.gca.linewidth,...
        'Units','inches',...
        'position',[public.plot.gca.ws_x  public.plot.gca.ws_y ...
     public.plot.gca.width public.plot.gca.height])
 
 %adjust the yticklabel
ax =gca;
ytick = ax.YTick;
expVal = -2;
new_label = cell(1,length(ytick));
for i =1:length(ytick)
    new_label{i}= sprintf('%g',ytick(i)/(10^expVal));
end
new_label{i}= '';
set(gca,'YTickLabel',new_label);
pos = get(gca,'Position');
 % offset = 0.00;  %how far to the right you want the exponent
 annotation('textbox',[0.0 0.9, 0.2, 0.2],... 
     'Units', 'pixels',...
    'String',['$\times 10^{' num2str(expVal) '}$'],...
    'fontsize',public.plot.gca.fontsize,...
    'fontname',public.plot.gca.fontname,...
    'Interpreter','latex',...
    'VerticalAlignment','bottom',...
    'EdgeColor','none')

%  pos = get(FIG1,'Position');
    set(FIG1, 'PaperUnits', 'inches','PaperPosition', ...
    [0  0 ...
    public.plot.gca.width  public.plot.gca.height],...
    'PaperSize', [public.plot.gca.width  public.plot.gca.height]);
name ='3D_BEB_case2_alg';
    picname1 = strcat(filepath,['\',name,'.pdf']);
% exportgraphics(FIG1,picname1,'ContentType','vector')
% saveas(FIG1,picname1)
%% two LCOs to plot

tspan = [0 50];
InitCond = LCO(:,1);
[tout,yout,yeout0,teout,yeout,ieout]=...
        Single_DS_Impacting_Hybrid_system_integration(A,R,C,InitCond,tspan,fs,equi_type);
    
FIG2 = figure;
set(FIG2,'units','inches')
h_ax2 = axes;

h2_1=plot3(yout(:,1)-1,yout(:,2),yout(:,3),'r--','linewidth',1.5,'displayname','Unstable LCO');
hold on; grid on
h2_0 = plot3(pe(1)-1,pe(2),pe(3),'ro','linewidth',2,'markerfacecolor','black','displayname','Stable PE');

% patch
box0_x = [0.034 0.049 0.049 0.034 0.034]*0;
box0_y = [-0.5 2.5 2.5 -0.5 -0.5 ];
box0_z = [-0.6 -0.6 0.6 0.6 -0.6 ]

V_0 = [box0_x(1), box0_y(1),box0_z(1);
    box0_x(2), box0_y(2),box0_z(2);
    box0_x(3), box0_y(3),box0_z(3);
    box0_x(4), box0_y(4),box0_z(4)];
F_0 = [1 2 3 4];


InitCond = LCO(:,2);
[tout,yout,yeout0,teout,yeout,ieout]=...
        Single_DS_Impacting_Hybrid_system_integration(A,R,C,InitCond,tspan,fs,equi_type);
 h2_2=plot3(yout(:,1)-1,yout(:,2),yout(:,3),'b-','linewidth',2,'displayname','Stable LCO');

% plot( T/pi,V_sign,'-','linewidth',1.2,'displayname','Sign of Velocity')

h2_3 =patch( 'Faces',F_0, 'Vertices',V_0,'EdgeColor',[0.1,0.1,0.3],'Linewidth',2, 'FaceColor', [0.1,0.1,0.3], 'FaceAlpha', 0.2,'displayname','$\Sigma$');
% xlim([0 5])
xlabel('$y_1$','interpreter','latex')
ylabel('$y_2$','interpreter','latex')
zlabel('$y_3$','interpreter','latex')

legend([h2_0 h2_1 h2_2 h2_3],'Location', 'best','interpreter','latex')

set(gca,'fontsize',public.plot.gca.fontsize,...
        'fontname',public.plot.gca.fontname,...
        'linewidth',public.plot.gca.linewidth,...
        'Units','inches',...
        'position',[public.plot.gca.ws_x  public.plot.gca.ws_y ...
     public.plot.gca.width public.plot.gca.height])
    set(FIG1, 'PaperUnits', 'inches','PaperPosition', ...
    [0  0 ...
    public.plot.gca.width  public.plot.gca.height],...
    'PaperSize', [public.plot.gca.width public.plot.gca.height]);
view(45,25)
%save 
    name ='3D_BEB_case2_LCO';
    picname2 = strcat(filepath,['\',name,'.pdf']);
% exportgraphics(FIG2,picname2,'ContentType','vector')



