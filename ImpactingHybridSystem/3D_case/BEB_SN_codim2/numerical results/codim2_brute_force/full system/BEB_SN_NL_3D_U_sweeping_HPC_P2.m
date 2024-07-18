%> ---------------- BRUTE FORCE DIAGRAM SIMULATION ----------------- %
clear
clc
close all

%% > initialization for the parallel computing
% addpath('/user/home/ib20968/Matlab_CodeStall/HPC/')
% run('HPC_add_path.m')
distcomp.feature( 'LocalUseMpiexec', true )

feature('numcores')

% local = parcluster('local')
% local.NumWorkers = 6;
% pool = local.parpool(6);

%%  define the parameter range interested --- 1D  parameter sweeping
T  = 12000;
fs = 512;  %> if fs = 0 will save the whole data 
t_kept = 200;
Brute_force_run_MIN_T   = 2000;
Brute_force_run_MAX_T   = T;
mfile_loc               = erase(mfilename('fullpath'),mfilename());
folder_2_save           = [mfile_loc,sprintf('BEB_SN_NL_3D_IHS_codim2_U_%s_P2',date)];
mkdir(folder_2_save);
if ~exist(folder_2_save, 'dir'); mkdir(folder_2_save); end

save_data_per_run = 0;

nordmark_vth            = 1e-6;
%
index 		=  [7;8]; %> the location of the flight velocity and the r
ds     		=  [1;0]; %> varying the parameter which controls the BEB
%
Par_ref = ...
        [-0.1,0.2,-0.5,1.78192697901049,1.6,4.79596124049338,0.025,0.0955685096979647]';
%> 
p_span          =[-0.05 0.025] ;
% >
p1_up_limit     = p_span(2);
p1_bot_limit    = p_span(1);
delta_1         = 0.0001;
p1_care_list    = p1_bot_limit:delta_1:p1_up_limit;


BEB_SN_NL_3D_U_sweeping;

copyfile([mfilename('fullpath'),'.m'],folder_2_save);
save(strcat(folder_2_save,['/',sprintf('BEB_SN_NL_3D_IHS_codim2_U_%s_P2',date)]))