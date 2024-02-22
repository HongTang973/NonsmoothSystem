%> plot the single case of the
close all

%%
load('brute_force_run_p_-0.0555_data.mat')
name1 ='3D_BEB_PD_codim2_portrait_0.0445.pdf';
figure
plot(tout, yout(:,1))
%
FIG1 = figure;
plot(yout(:,1), yout(:,2),'k-','LineWidth',1.0)
hold on
xlabel('$y_1$','Interpreter','latex')
ylabel('$y_2$','Interpreter','latex')
% xlim([-0.025 0.31])
adj_plot_theme_II(FIG1)
grid on
ax =gca;
ax.GridLineWidth = 1;
ytick = ax.YTick;
% expVal = -3;
new_label = cell(1,length(ytick));
for i =2:2:length(ytick)
    new_label{i}= sprintf('%g',ytick(i));
end
set(gca,'YTickLabel',new_label);
%>
xtick = ax.XTick;
x_new_label = cell(1,length(xtick));
for i =1:2:length(xtick)
    x_new_label{i}= sprintf('%g',xtick(i));
end
set(gca,'XTickLabel',x_new_label);

% exportgraphics(FIG1, name1,'ContentType','vector')
%%
name2 ='3D_BEB_PD_codim2_portrait_0.028.pdf';
load('brute_force_run_p_-0.0972_data.mat')
FIG2 = figure;
plot(yout(:,1), yout(:,2),'k-','LineWidth',1.0)
hold on
xlabel('$y_1$','Interpreter','latex')
ylabel('$y_2$','Interpreter','latex')
% xlim([-0.025 0.31])
adj_plot_theme_II(FIG2)
grid on
ax =gca;
ax.GridLineWidth = 1;
ytick = ax.YTick;
expVal = -3;
new_label = cell(1,length(ytick));
for i =2:2:length(ytick)
    new_label{i}= sprintf('%g',ytick(i)*10^-expVal);
end
set(gca,'YTickLabel',new_label);
%>
xtick = ax.XTick;
expVal = -3;
x_new_label = cell(1,length(xtick));
for i =1:2:length(xtick)
    x_new_label{i}= sprintf('%g',xtick(i)*10^-expVal);
end
set(gca,'XTickLabel',x_new_label);

annotation('textbox',[0.00 0.85, 0.2, 0.2],...
    'Units', 'pixels',...
    'String',['$\times 10^{' num2str(expVal) '}$'],...
    'fontsize',16,...
    'fontname','times new roman',...
    'Interpreter','latex',...
    'VerticalAlignment','bottom',...
    'EdgeColor','none')
%
annotation('textbox',[0.80 0.00, 0.2, 0.2],...
    'Units', 'pixels',...
    'String',['$\times 10^{' num2str(expVal) '}$'],...
    'fontsize',16,...
    'fontname','times new roman',...
    'Interpreter','latex',...
    'VerticalAlignment','bottom',...
    'EdgeColor','none')

% exportgraphics(FIG2, name2,'ContentType','vector')


%% 
name3 ='3D_BEB_PD_codim2_portrait_0.027.pdf';
load('brute_force_run_p_-0.0973_data.mat')
FIG3 = figure;
plot(yout(:,1), yout(:,2),'k-','LineWidth',1.0)
hold on
xlabel('$y_1$','Interpreter','latex')
ylabel('$y_2$','Interpreter','latex')
xlim([-1 20]*1e-3)
adj_plot_theme_II(FIG3)
grid on
ax =gca;
ax.GridLineWidth = 1;
ytick = ax.YTick;
expVal = -3;
new_label = cell(1,length(ytick));
for i =2:2:length(ytick)
    new_label{i}= sprintf('%g',ytick(i)*10^-expVal);
end
set(gca,'YTickLabel',new_label);
%>
xtick = ax.XTick;
expVal = -3;
x_new_label = cell(1,length(xtick));
for i =1:2:length(xtick)
    x_new_label{i}= sprintf('%g',xtick(i)*10^-expVal);
end
set(gca,'XTickLabel',x_new_label);

annotation('textbox',[0.00 0.85, 0.2, 0.2],...
    'Units', 'pixels',...
    'String',['$\times 10^{' num2str(expVal) '}$'],...
    'fontsize',16,...
    'fontname','times new roman',...
    'Interpreter','latex',...
    'VerticalAlignment','bottom',...
    'EdgeColor','none')
%
annotation('textbox',[0.80 0.00, 0.2, 0.2],...
    'Units', 'pixels',...
    'String',['$\times 10^{' num2str(expVal) '}$'],...
    'fontsize',16,...
    'fontname','times new roman',...
    'Interpreter','latex',...
    'VerticalAlignment','bottom',...
    'EdgeColor','none')

% exportgraphics(FIG3, name3,'ContentType','vector')


%% 

name4 ='3D_BEB_PD_codim2_portrait_0.015.pdf';
load('brute_force_run_p_-0.0985_data.mat')
FIG4 = figure;
plot(yout(:,1), yout(:,2),'k-','LineWidth',1.0)
hold on
xlabel('$y_1$','Interpreter','latex')
ylabel('$y_2$','Interpreter','latex')
xlim([-1 11]*1e-3)
adj_plot_theme_II(FIG4)
grid on
ax =gca;
ax.GridLineWidth = 1;
ytick = ax.YTick;
expVal = -3;
new_label = cell(1,length(ytick));
for i =2:2:length(ytick)
    new_label{i}= sprintf('%g',ytick(i)*10^-expVal);
end
set(gca,'YTickLabel',new_label);
%>
xtick = ax.XTick;
expVal = -3;
x_new_label = cell(1,length(xtick));
for i =1:2:length(xtick)
    x_new_label{i}= sprintf('%g',xtick(i)*10^-expVal);
end
set(gca,'XTickLabel',x_new_label);

annotation('textbox',[0.00 0.85, 0.2, 0.2],...
    'Units', 'pixels',...
    'String',['$\times 10^{' num2str(expVal) '}$'],...
    'fontsize',16,...
    'fontname','times new roman',...
    'Interpreter','latex',...
    'VerticalAlignment','bottom',...
    'EdgeColor','none')
%
annotation('textbox',[0.80 0.00, 0.2, 0.2],...
    'Units', 'pixels',...
    'String',['$\times 10^{' num2str(expVal) '}$'],...
    'fontsize',16,...
    'fontname','times new roman',...
    'Interpreter','latex',...
    'VerticalAlignment','bottom',...
    'EdgeColor','none')

annotation('textarrow',[0.24 0.30],[0.62 0.62])
rectangle('Position',[-2e-4 -1e-3  4e-4 2e-3],'LineWidth',1,'EdgeColor',[1 0 0])

axes('position',[.30 .55 .15 .15])
% box on % put box around new pair of axes
plot(yout(:,1), yout(:,2),'k-','LineWidth',1.0)
grid on
xlim([-1e-4 1e-4])
ylim([-1e-3 1e-3])
set(gca, 'color', 'none', 'LineWidth',1,'Yticklabel',[],'FontName','Times New Roman');
% set(gca,'Xticklabel',[])

% exportgraphics(FIG4, name4,'ContentType','vector')
