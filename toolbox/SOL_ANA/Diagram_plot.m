%% post processing and plots
close all
% chose bifurcation parameter
para=menu('Which parameter was used for analysis?','Velocity','Damping',...
            'Preload','Gapsize','AoA');
      
        switch para
            case 1    
                if exist('par_V_list','var')
                    V_care_list=par_V_list;
                   
                elseif ~exist('V_care_list','var')
                    
                V_care_list=10:0.1:90;
                end
               bifur_para=V_care_list; 
               
               
%                xlabel_txt='Velocity $U/\omega_ab$';
%                non_dimension_factor=1;
%                
                xlabel_txt='$\epsilon$';
               non_dimension_factor=1;
               
            case 2
               xlabel_txt='Damping $D/D_0$';
               non_dimension_factor=1;
            case 3
               xlabel_txt='Preload $\times 1 rad \cdot k_{beta}$';
               non_dimension_factor=1;
            case 4
               xlabel_txt='Gap size $\delta/0.01$';
               non_dimension_factor=100;
            case 5
               xlabel_txt='Angle of attack \, / rad';
               non_dimension_factor=1;
        end
% baoyuanxiang 
% beta
n_col=max(beta_length_data);
beta_diagram=zeros(length(beta_length_data),n_col);
temp_indx=1;
for j=1:length(beta_length_data)
    beta_diagram(j,1:beta_length_data(j))=beta_diogram_matrix(temp_indx:temp_indx+beta_length_data(j)-1);
    beta_diagram(j,beta_length_data(j)+1:end)=nan;
    temp_indx=temp_indx+beta_length_data(j);
end
figure(1)
hold on
% Peter_CreatePlotInOrigin(beta_diagram,'beta_diagram')   
plot(beta_diagram(:,1)*non_dimension_factor,beta_diagram(:,2:end),'k.','markersize',4)
hold on


xlabel(xlabel_txt,'interpreter','latex');
ylabel('Flap motion amplitude$/\delta$','interpreter','latex');

% xlim([0.6 0.74])
%  ylim([-1.6 1.5])
title('Flap motion ');
set(gca,'fontname','Times New Roman','fontsize',12,'linewidth',1.2)

% AOA
n_col=max(alpha_length_data);
alpha_diagram=zeros(length(alpha_length_data),n_col);
temp_indx=1;
for j=1:length(alpha_length_data)
    alpha_diagram(j,1:alpha_length_data(j))=alpha_diogram_matrix(temp_indx:temp_indx+alpha_length_data(j)-1);
    alpha_diagram(j,alpha_length_data(j)+1:end)=nan;
    temp_indx=temp_indx+alpha_length_data(j);
end
figure(2)
plot(alpha_diagram(:,1)*non_dimension_factor,alpha_diagram(:,2:end),'k.','markersize',4)
hold on

xlabel(xlabel_txt,'interpreter','latex');
ylabel('Pitch motion amplitude/rad');

% xlim([0.6 0.74])
% ylim([0.026 -0.964])
% ylim([])
title('Pitch motion ');
set(gca,'fontname','Times New Roman','fontsize',12,'linewidth',1.2)

% Heave
n_col=max(heave_length_data);
heave_diagram=zeros(length(heave_length_data),n_col);
temp_indx=1;
for j=1:length(heave_length_data)
    heave_diagram(j,1:heave_length_data(j))=heave_diogram_matrix(temp_indx:temp_indx+heave_length_data(j)-1);
    heave_diagram(j,heave_length_data(j)+1:end)=nan;
    temp_indx=temp_indx+heave_length_data(j);
end

% 
figure(3)
plot(heave_diagram(:,1)*non_dimension_factor,heave_diagram(:,2:end),'k.','markersize',4)
hold on

xlabel(xlabel_txt,'interpreter','latex');
ylabel('Nondimensional heave motion h/b');

% xlim([0.6 0.74])
% ylim([-0.025 0.005])
title('Heave motion ');
set(gca,'fontname','Times New Roman','fontsize',12,'linewidth',1.2)
%
%

% 

%
%
%% adding
% figure(1)
% h1=plot(beta_diagram(:,1)*non_dimension_factor,beta_equi_matrix,'--','linewidth',1);legend(h1,'Static equilibrium')
% figure(2)
% h2=plot(beta_diagram(:,1)*non_dimension_factor,alpha_equi_matrix,'--','linewidth',1);legend(h2,'Static equilibrium')
% figure(3)
% h3=plot(beta_diagram(:,1)*non_dimension_factor,heave_equi_matrix,'--','linewidth',1);legend(h3,'Static equilibrium')
