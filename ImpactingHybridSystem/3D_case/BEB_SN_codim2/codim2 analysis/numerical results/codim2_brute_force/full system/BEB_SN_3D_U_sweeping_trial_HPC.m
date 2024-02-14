
%> ---------------- BRUTE FORCE DIAGRAM SIMULATION ----------------- %
clear
clc
close all

%% > initialization for the parallel computing
% addpath('/user/home/ib20968/Matlab_CodeStall/HPC/')
% run('HPC_add_path.m')
Initialization_for_Par;

%%  define the parameter range interested --- 1D  parameter sweeping
T  = 80;
fs = 512;  %> if fs = 0 will save the whole data 
t_kept = 5;
Brute_force_run_MIN_T = 50;
Brute_force_run_MAX_T = T;
mfile_loc =erase(mfilename('fullpath'),mfilename());
folder_2_save = [mfile_loc,sprintf('Airfoil_model_PD_codim2_U_%s_P1',date)];
if ~exist(folder_2_save, 'dir'); mkdir(folder_2_save); end
save_data_per_run = 1;
%
index 		=  [1;5]; %> the location of the flight velocity and the damping
ds     		=  [1;0]; %> varying the parameter which controls the BEB
%
Par_ref = ...
        [19.4300;0.62;1;1;1;0.133314113099967];
    
%> 
V_range 	= [19.2 19.8];
p_span      = V_range - Par_ref(1);
% >
p1_up_limit     = p_span(2);
p1_bot_limit    = p_span(1);
delta_1         = 0.2;
p1_care_list    = p1_bot_limit:delta_1:p1_up_limit;
%>

BEB_SN_3D_U_sweeping;

copyfile([mfilename('fullpath'),'.m'],folder_2_save);
save(strcat(folder_2_save,['\',sprintf('Airfoil_model_SN_codim2_U_%s_P1.mat',date)]))