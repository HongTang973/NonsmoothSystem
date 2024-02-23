%> plot the single case of the
close all

%%
% load('brute_force_run_p_-0.0555_data.mat')
% name1 ='3D_BEB_PD_codim2_portrait_0.0445.pdf';
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
