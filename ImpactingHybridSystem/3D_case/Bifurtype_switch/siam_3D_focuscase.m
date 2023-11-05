clc
close all
clear
% 3D case numerical playing around
% define the public contaning the common features
public.plot.addmissible_ob.color = [0 ,0 , 0]; %black
public.plot.addmissible_ob.linewidth = 1;
public.plot.addmissible_ob.style= '-';
%
public.plot.S_po.color = [0 ,0 , 1]; % blue
public.plot.S_po.linewidth = 2;
public.plot.S_po.style = '-';

%
public.plot.US_po.linewidth = 2;
public.plot.US_po.color = [1,0,0];
public.plot.US_po.style = '--';
%

public.plot.gca.linewidth = 3;
public.plot.gca.fontname = 'Times New Roman';
public.plot.gca.fontsize = 16;
public.plot.gca.ws_x = 0.9;
public.plot.gca.ws_y = 0.7;
public.plot.gca.width = 4.5208;
public.plot.gca.height = 3.5656;
% if equi_type == 'admissible'
%         pre = -1;
%  elseif equi_type == 'pseudo'
%         pre = 1;
%  end
equi_type = 1;

filepath = [pwd,'\figures'];
if ~exist(filepath, 'dir')
       mkdir(filepath)
end
%% use nondimensional egenvalues we can get to two cases with
% lambda_3 = 1 or lambda_3 = -1
% firstly with the -1 case

% defined in the OB form with special case
alpha = -10;
beta =10;

p = struct();
p.lambda_1 = alpha + beta*1j;
p.lambda_2 =  alpha - beta*1j;
p.lambda_3 = 1;
p.b1 = 0;
p.b2 = 2.5;
p.b3 = 1;
p.phi = pi/2;
p.theta = 0;
% then transfer to general Jordan form
[A,B,C,phi2,theta2] = OB2Jordan(p)

%
% get from OB2Jordan
p01 = struct();
p01.lambda_1 =  alpha + beta*1j;
p01.lambda_2 =  alpha - beta*1j;
p01.lambda_3 = 1;
p01.phi = phi2;
p01.theta = theta2;

p01.b1 = B(1);
p01.b2 = B(2);
p01.b3 = B(3);

% The second case
p02 = p01;
p02.lambda_3 = -1;

% transform to the algorithm compatible form
% case 1 -- node case
[A1,B1,C1,R1,A_s1]= Matrices_3D_IHS_Gform(p01,'Case2')
[A2,B2,C2,R2,A_s2]= Matrices_3D_IHS_Gform(p02,'Case2')
as=eig(A_s1)
%
[V1,D1]=eig(A1);

% determine the sampling frequency
Omega = max(abs(diag(D1)));
fs = 200*ceil(2*(Omega/2/pi));
% for different Evolution time
a =0;
b = 0.04;
delta = 0.001;
T = a*pi:delta:b*pi;

%
MAX1 = zeros(1,length(T));
F_1 = zeros(1,length(T));
V_sign1 = zeros(1,length(T));
LOCI1= zeros(size(A,1),length(T));
vector1= zeros(size(A,1),length(T));
%
MAX2 = zeros(1,length(T));
F_2 = zeros(1,length(T));
V_sign2 = zeros(1,length(T));
LOCI2= zeros(size(A,1),length(T));
vector2= zeros(size(A,1),length(T));

for i=1:length(T)
    [V_sign1(i),LOCI1(:,i),MAX1(i),vector1(:,i),F_1(i),~] = LCO_Det_search(T(i),R1,A1,C1,equi_type);
     [V_sign2(i),LOCI2(:,i),MAX2(i),vector2(:,i),F_2(i),~] = LCO_Det_search(T(i),R2,A2,C2,equi_type);
end

%% make the plot


FIG1 = figure;set(FIG1,'Units','inches');
hax1 = axes;
h11 = plot(T/pi,F_1,'r-','linewidth',1.2,'displayname','$Case II - (i)$');
hold on
h13 = plot(T/pi,F_2,'k-','linewidth',1.2,'displayname','$Case II - (ii)$');
%  plot( T/pi,V_sign,'-','linewidth',1.2,'displayname','Sign of Velocity')
% h10 = plot([a b],[0 0],'b-','displayname','zero line');
% legend

xlabel('\tau /\pi')
ylabel('$p(t)$','interpreter','latex')
xlim([0 b])

grid on
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
    new_label{i}= sprintf('%.1g',ytick(i)/(10^expVal));
end
set(gca,'YTickLabel',new_label);

 offset = 0.00;  %how far to the right you want the exponent
 annotation('textbox',[0.00+offset 0.9, 0.2, 0.2],... 
     'Units', 'pixels',...
    'String',['$\times 10^{' num2str(expVal) '}$'],...
    'fontsize',public.plot.gca.fontsize,...
    'fontname',public.plot.gca.fontname,...
    'Interpreter','latex',...
    'VerticalAlignment','bottom',...
    'EdgeColor','none')
%save
legend([ h11,  h13],'location','best','interpreter','latex')
%
pos = get(FIG1,'Position');
set(FIG1, 'PaperUnits', 'inches','PaperPosition', ...
    [0  0 ...
    public.plot.gca.width  public.plot.gca.height],...
    'PaperSize', pos(3:4));
%
name ='3DOF_switchcase2';
    picname1 = strcat(filepath,['\',name,'.pdf']);
   exportgraphics(FIG1,picname1,'ContentType','vector')