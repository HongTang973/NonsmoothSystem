
% Codimiension one continuation parameter analysis
%  use continuation to track the L vs q
clc
close all
clear

%
rootpath = erase(mfilename('fullpath'),mfilename());
filepath = [rootpath,'\figures'];
if ~exist(filepath, 'dir')
       mkdir(filepath)
end
%

public.plot.gca.linewidth = 3;
public.plot.gca.fontname = 'Times New Roman';
public.plot.gca.fontsize = 16;
public.plot.gca.ws_x = 0.9;
public.plot.gca.ws_y = 0.7;
public.plot.gca.width = 4.5208;
public.plot.gca.height = 3.5656;


%% in order to make the  two ways of passing varible compatible in 
% this process
% first; we give a dict  par_dict containing keys and values 
% this par_dict can be converted to structure via 
%  p = entries(par_dict, "struct")
% and can be converted to a table by 
% table  = entries(par_dict);




par = [343.14, 0.78, 0.05 , 1.058,  100, 1]';
keys = {'a', 'rx', 'kappa','L', 'RT','xg'};
% record the index 
S.index = 1:length(par);

par_struct = struct_construction(keys, par);
p               = nd_QWM_N1_gas(par_struct);
% clarify the equilibrium and the re-construnction order

x_list          = [0:0.001:2]*p.x_max/p.x_ref;
l_cr            = p.L_cr(x_list);  


ind_p  =    6; ind_x  =    4;
p_span      = [0 1];
x_span      = [0 2];

prob.par    =    par;
prob.keys   =    keys;
prob.ind_p  =    6;
prob.ind_x  =    4;
prob.dir    =    1;
prob.p_span = [0 1];
prob.x_span = [0 2];
prob.ZeroFunctions = @ZeroFunctions_Valve_QWM_HB;

%

% 

tic
prob_1 = prob;
prob_1.p_span = [1e-2 1];
prob_1.co_dh  = 1e-2;
[u1,iter1] = codim1_PC(prob_1);
% [u1,iter1] = codim1_PC(@ZeroFunctions_Valve_QWM_HB,par,ind_p,ind_x,p_span,x_span,1);
    
prob_2 = prob;
prob_2.par    = [343.140000000000;0.780000000000000;0.0500000000000000;0.0756309428606699;100;0.00198124891183716];
% prob_2.par=[343.140000000000;0.780000000000000;0.0500000000000000;0.0896181209608336;100;0.00446422825123808];
% prob_2.par      = u1(:,end);
prob_2.p_span = [1e-4 1e-2];
prob_2.co_dh  = 1e-4;
[u2,iter2] = codim1_PC(prob_2);
% 
prob_3          = prob;
% prob_3.par      = [343.140000000000;0.780000000000000;0.0500000000000000;0.0756309428606699;100;0.00198124891183716];
prob_3.par      = u2(:,end);
prob_3.p_span   = [1e-6 1e-4];
prob_3.co_dh    = 1e-6;
[u3,iter3]      = codim1_PC(prob_3);
%
prob_4          = prob;
prob_4.par      = u3(:,end);
prob_4.p_span   = [1e-8 1e-6];
prob_4.co_dh    = 1e-8;
[u4,iter4]      = codim1_PC(prob_4);
% % from the liquid side
% par = [343.14, 0.78, 0.05 , 4.003,  100, 0]';
prob_5          = prob;
prob_5.par      = [343.14, 0.78, 0.05 , 4.003,  100, 0]';
prob_5.p_span   = [0 1e-6];
prob_5.x_span = [0 5];
prob_5.dir      = 1;
prob_5.co_dh    = 1e-10;
[u5,iter5]      = codim1_PC(prob_5);
% 

% [u3,iter3] = codim1_PC(@ZeroFunctions_Valve_QWM_HB,par,ind_p,ind_x,p_span,x_span,-1);
CPU_tim = toc
figure
plot(u1(prob.ind_p,:),u1(prob.ind_x,:))
hold on
% figure
plot(u2(prob.ind_p,:),u2(prob.ind_x,:))

FIG1 = figure;set(FIG1,'Units','inches');
hax1 = axes;
h10 = plot(u1(ind_p,:),u1(ind_x,:),'color','black','linewidth',1.5,'displayname','CO');
hold on
plot(u2(ind_p,:),u2(ind_x,:),'color','blue','linewidth',1.5)
plot(u1(ind_p,1),u1(ind_x,1),'o','color','red','linewidth',2,'markersize',8,'markerfacecolor',[1 1 1])
% 
h20 =plot(x_list*p.x_ref/p.x_max*100, l_cr,'r--', 'linewidth',1.5, 'displayname','QWM-Ana');
legend([h10, h20], 'location', 'best')
xlabel('$xg$','interpreter','latex')
ylabel('$L_{cr}$','interpreter','latex')

set(gca,'fontsize',public.plot.gca.fontsize,...
        'fontname',public.plot.gca.fontname,...
        'linewidth',public.plot.gca.linewidth,...
        'Units','inches',...
        'position',[public.plot.gca.ws_x  public.plot.gca.ws_y ...
     public.plot.gca.width public.plot.gca.height])
 pos = get(FIG1,'Position');
    set(FIG1, 'PaperUnits', 'inches','PaperPosition', ...
    [0  0 ...
    public.plot.gca.width  public.plot.gca.height],...
    'PaperSize', pos(3:4));
 %save 
%     name ='case1_continuation_r_T';
%     picname1 = strcat(filepath,['\',name,'.pdf']);
% %    exportgraphics(FIG1,picname1,'ContentType','vector')
%    % save the pics 
%     saveas(FIG1,picname1)
%%  monitoring functions to track the cared feature value
% monitor = [Mono_p;Salt_p;D_s;ve;sub_norm];

len_1 = size(u1,2);
len_2 = size(u2,2);
Obs_1 = [];
Obs_2 = [];
for i=1:len_1
    par_temp = u1(:,i);
    [monitor] = Obeservation_Funs_Valve_QWM(par_temp);
    Obs_1 = [Obs_1,monitor];
end

for j=1:len_2
    par_temp = u2(:,j);
    [monitor] = Obeservation_Funs_Valve_QWM(par_temp);
    Obs_2 = [Obs_2,monitor];
end

% track the greatest floquet multiplier
delta       =  0.01;
gfm         = [flip(Obs_1(10,2:end)),Obs_2(10,:)];

p_list      = [flip(u1(ind_p,2:end)),u2(ind_p,:)];
x_list      = [flip(u1(ind_x,2:end)),u2(ind_x,:)];

disc2unit   = (gfm -1);
% there is a/some cross points
index0 =sign(disc2unit);
index1 = abs(diff(index0))>0;
% filter the singularity case
index_s = abs(diff(disc2unit))/delta < (1/delta);
index1 = index1 & index_s;
index2 = [index1,0];
index3 = [0,index1];
index2 = find(index2==1);
index3 = find(index3==1);
% section partia
ratio       = abs(disc2unit(index2))./(abs(disc2unit(index2))+abs(disc2unit(index3)));
p_chosen    = (1-ratio).*p_list(index2)+ratio.*p_list(index3);
MAX_chosen  = (1-ratio).*gfm(index2)+ratio.*gfm(index3);
F1_chosen   = (1-ratio).*disc2unit(index2)+ratio.*disc2unit(index3);
% 
seg_1_p = p_list(1:index3(1));
seg_2_p = p_list(index3(1):index3(2));
seg_3_p = p_list(index3(2):index3(3));
seg_4_p = p_list(index3(3):end);
seg_1_gfm = gfm(1:index3(1));
seg_2_gfm = gfm(index3(1):index3(2));
seg_3_gfm = gfm(index3(2):index3(3));
seg_4_gfm = gfm(index3(3):end);

seg_1_x = x_list(1:index3(1));
seg_2_x = x_list(index3(1):index3(2));
seg_3_x = x_list(index3(2):index3(3));
seg_4_x = x_list(index3(3):end);

figure;
plot(seg_1_p,seg_1_x, 'r')
hold on
plot(seg_2_p,seg_2_x, 'g')
plot(seg_3_p,seg_3_x, 'k')
plot(seg_4_p,seg_4_x, 'b')

FIG2 = figure;
hax2 =axes;
plot(u1(ind_p,:),abs(Obs_1(10,:)),'.','color','black')
hold on
plot(u2(ind_p,:),abs(Obs_2(10,:)),'.','color','blue')

[bif_val,bif_ind] = min(u1(ind_p,:));

plot(hax1,bif_val,u1(end,bif_ind),'o','color','green','linewidth',2,'markersize',8,'markerfacecolor',[1 1 1])
plot(hax1,p_list(index3), x_list(index3),'o','color','blue','linewidth',2,'markersize',8,'markerfacecolor',[1 1 1])

%% unstable branch and stable branch

u1_s = u1(:,1:bif_ind);

u1_us = u1(:,bif_ind:end);
FIG3 = figure;set(FIG3,'Units','inches');

h3_1=plot(seg_2_p,seg_2_x,'color','blue','linewidth',2,'displayname','Stable branch');
hold on
plot(seg_4_p,seg_4_x,'color','blue','linewidth',2)
h3_2=plot(seg_1_p,seg_1_x,'--','color','red','linewidth',2,'displayname','Unstable branch');
h3_0 = plot(u1(ind_p,1),u1(end,1),'o','color','black','linewidth',2,'markersize',8,'markerfacecolor',[1 1 1],'displayname',['$\mbox{Current point:}$' newline '$r=',sprintf('%3g',par(2)),'$']);
h3_4 = plot(bif_val,u1(end,bif_ind),'o','color','green','linewidth',2,'markersize',8,'markerfacecolor',[1 1 1],'displayname',['$\mbox{Turning point:}$' newline '$r \approx',sprintf('%3g',bif_val),'$']);
% '$\vbox{\noindent \mbox{Current point: }\newline r=0.72}$'
h3_3 = plot(seg_2_p(1),seg_2_x(1),'o','color','red','linewidth',2,'markersize',8,'markerfacecolor',[1 1 1],'displayname',['$\mbox{B point:}$' newline '$r \approx',sprintf('%3g',seg_2_p(1)),'$']);
xlabel('$r$','interpreter','latex')
ylabel('$T^*$','interpreter','latex')
ylim([-0.01 0.7])

legend([h3_0 h3_1 h3_3 h3_4 h3_2],'location','best','interpreter','latex')
set(gca,'fontsize',public.plot.gca.fontsize,...
        'fontname',public.plot.gca.fontname,...
        'linewidth',public.plot.gca.linewidth,...
        'Units','inches',...
        'position',[public.plot.gca.ws_x  public.plot.gca.ws_y ...
     public.plot.gca.width public.plot.gca.height])
 pos = get(FIG3,'Position');
    set(FIG3, 'PaperUnits', 'inches','PaperPosition', ...
    [0  0 ...
    public.plot.gca.width  public.plot.gca.height],...
    'PaperSize', pos(3:4));
 %save 
    name ='case1_co_r_T_US';
    picname1 = strcat(filepath,['\',name,'.pdf']);
%    exportgraphics(FIG3,picname1,'ContentType','vector')
   % save the pics 
    saveas(FIG3,picname1)





