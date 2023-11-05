clc
clear
close all
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
p = struct();
p.lambda_1 = -1.2;
p.lambda_2 = -0.01;
p.lambda_3 = -1;
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
p01.lambda_1 = -1.2;
p01.lambda_2 = 0.01;
p01.lambda_3 = -1;
p01.phi = phi2;
p01.theta = theta2;

p01.b1 = B(1);
p01.b2 = B(2);
p01.b3 = B(3);

% The second case
p02 = p01;
p02.lambda_2 = -0.01;

% transform to the algorithm compatible form
% case 1 -- node case
[A1,B1,C1,R1,A_s1]= Matrices_3D_IHS_Gform(p01,'Case1')
[A2,B2,C2,R2,A_s2]= Matrices_3D_IHS_Gform(p02,'Case1')
as=eig(A_s1)
%
[V1,D1]=eig(A1);

% determine the sampling frequency
Omega = max(abs(diag(D1)));
fs = 200*ceil(2*(Omega/2/pi));
% for different Evolution time
a =0;
b = 1;
delta = 0.005;
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
    [V_sign1(i),LOCI1(:,i),MAX1(i),vector1(:,i),F_1(i)] = LCO_Det_search(T(i),R1,A1,C1,equi_type);
     [V_sign2(i),LOCI2(:,i),MAX2(i),vector2(:,i),F_2(i)] = LCO_Det_search(T(i),R2,A2,C2,equi_type);
end

%% make the plot


FIG1 = figure;set(FIG1,'Units','inches');
hax1 = axes;
h11 = plot(T/pi,F_1,'r-','linewidth',1.2,'displayname','$Case I - (i)$');
hold on
h13 = plot(T/pi,F_2,'k-','linewidth',1.2,'displayname','$Case I - (ii)$');
%  plot( T/pi,V_sign,'-','linewidth',1.2,'displayname','Sign of Velocity')
% h10 = plot([a b],[0 0],'b-','displayname','zero line');
% legend

xlabel('\tau /\pi')
ylabel('$p(t)$','interpreter','latex')
xlim([0 0.4])

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
expVal = -3;
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
name ='3DOF_switchcase1';
    picname1 = strcat(filepath,['\',name,'.pdf']);
   exportgraphics(FIG1,picname1,'ContentType','vector')
%    saveas(FIG1,picname1)
index0 =sign(F_1);
index1 = abs(diff(index0))>0;
% filter the singularity case
index_s = abs(diff(F_1))/delta < (1/delta);
index1 =index1 & index_s;
index2 = [index1,0];
index3 = [0,index1];
index2=find(index2==1);
index3=find(index3==1);
% section partia
ratio = abs(F_1(index2))./(abs(F_1(index2))+abs(F_1(index3)));
T_chosen =(1-ratio).*T(index2)+ratio.*T(index3);
MAX_chosen = (1-ratio).*MAX(index2)+ratio.*MAX(index3);
F1_chosen = (1-ratio).*F_1(index2)+ratio.*F_1(index3);
Sign_chosen =(1-ratio).*V_sign(index2)+ratio.*V_sign(index3);
%          figure(1)
% plot(T_chosen/pi,Sign_chosen,'bo')
%     %
index4 =  (abs(F1_chosen)<1e-3) & (T_chosen>1e-3);
% index4 = (Sign_chosen == 1) & (abs(F1_chosen)<1e-3) %& (abs(MAX_chosen)<1e-3);
T_chosen  =T_chosen(index4);
F1_chosen =F1_chosen(index4);

% we know there is a LCO

% get the LCO with statevariables
[V_sign_LCO1,~,F11,LCO1,~] = LCO_Det_search(T_chosen,R,A,C,equi_type);

h12 = plot(hax1,T_chosen/pi,F1_chosen,'bo','displayname',['the ',num2str(1),'st candidate']);

legend([h10, h11, h12],'location','best','interpreter','latex')

V_sign_LCO


%check the floque multipliers of the foud LCO
[Mono_p,Salt_p]=Floque_Multipliers(T_chosen,LCO1,R1,A1,C1)
if max(abs(Salt_p))>1+1e-6
    disp([' LCO1  is unstable!'])
else
    disp([' LCO1  is stable!'])
end

% for the second case
%check the floque multipliers of the foud LCO
[Mono_p,Salt_p]=Floque_Multipliers(T_chosen,LCO2,R1,A1,C1)
if max(abs(Salt_p))>1+1e-6
    disp([' LCO1  is unstable!'])
else
    disp([' LCO1  is stable!'])
end

keyboard
%
tspan = [0 50];
if ~isempty(T_chosen)
    
    for i =1:length(T_chosen)
        InitCond = LCO(:,i);
        [tout,yout,yeout0,teout,yeout,ieout]=...
            Single_DS_Impacting_Hybrid_system_integration(A,R,C,V_sign_LCO(i)*InitCond,tspan,fs,V_sign_LCO(i)*equi_type);
        eval(['tout',num2str(i),'=tout;']);
        eval(['yout',num2str(i),'=yout;']);
        eval(['yeout0',num2str(i),'=yeout0;']);
        eval(['teout',num2str(i),'=teout;']);
        eval(['yeout',num2str(i),'=yeout;']);
        eval(['ieout',num2str(i),'=ieout;']);
        eval(['delta_t',num2str(i),'=teout',num2str(i),'(2:end)-teout',num2str(i),'(1:end-1);']);
        %Period
        eval(['Period',num2str(i),'=delta_t',num2str(i),'/pi;']);
        figure
        plot(tout/pi,yout(:,C>0),'r-','linewidth',1.4)
        title(['LCO ',num2str(i),'''s stability'])
        set(gca,'fontname','times new roman','fontsize',12)
        xlabel('t/s')
        figure
        plot3(yout(:,1),yout(:,2),yout(:,3),'k-','linewidth',1.2)
        title(['LCO ',num2str(i),'''s phase portrait'])
        set(gca,'fontname','times new roman','fontsize',12)
        grid on
        figure
        plot(yout(:,1),yout(:,2),'k-','linewidth',1.2)
        title(['LCO ',num2str(i),'''s phase portrait'])
        set(gca,'fontname','times new roman','fontsize',12)
        grid on
        figure
        plot(yout(:,2),yout(:,3),'k-','linewidth',1.2)
        title(['LCO ',num2str(i),'''s phase portrait'])
        set(gca,'fontname','times new roman','fontsize',12)
        grid on
    end
end

[tout,yout,yeout0,teout,yeout,ieout]=...
    Single_DS_Impacting_Hybrid_system_integration(A,R,C,[1;23;-300],[0 1],fs,equi_type);figure;plot(tout,yout(:,1))

