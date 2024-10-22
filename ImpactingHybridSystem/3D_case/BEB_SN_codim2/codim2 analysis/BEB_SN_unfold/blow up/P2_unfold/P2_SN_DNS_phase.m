%> simulate the time history at the SN point.
%> ---------------- BRUTE FORCE DIAGRAM SIMULATION ----------------- %
t_kept =100;
Brute_force_run_MIN_T =1000;
Brute_force_run_MAX_T = 5000;
fs                    =32;
save_data_per_run = 0;
equi_type         = 1;
point_1 = [-0.1,0.2,-0.525,1.78192697901049+0.0955685096979647,1.6,4.79596124049338]';
p_care  = 0;
index   = [3,4];
ds      = [1;0];
%> initialize the
A_diogram_matrix=[]; % vector to keep result
A_length_data=[];    % record length for per result in order to recover
% data from vector to matrix
%
B_diogram_matrix=[];
B_length_data=[];
%
run_info_collection     =[];
%
C_diogram_matrix=[];
C_length_data=[];

%> start the computation using paralleledcomputation
tic


disp(p_care)
%> choose one point with fixed eta
point_1(index) = point_1(index) + p_care*ds;
%
[A,B,C,R, T_2_det]  = par2NForm_Lienard(point_1);
% x_00                = IC_generaotr(point_1(6),R,A,C,equi_type)
prob = IHS_3D_par2prob(point_1,equi_type);
%

[tout,yout,yeout0,teout,yeout,ieout,run_info]=...
    brute_force_blow_up_Shilnikov_DNS_SD_int_solver(prob,t_kept,Brute_force_run_MIN_T, Brute_force_run_MAX_T);

%> initialize the buffer variables for single run
%% >  save the data per run
if save_data_per_run
    %save a small file with the bifurcation variables for this run
    file_name  = sprintf('brute_force_run_p_%g_data.mat',p_care);
    mat_name   = strcat(folder_2_save,['//',file_name]);
    vnames     ={  'tout','yout','yeout0','teout','yeout','ieout', 'run_info'};
    parsave_named(mat_name,vnames,tout,yout,yeout0,teout,yeout,ieout,run_info);
    %  note-- parsave_named is needed because the usual matlab save
    %  function doesn't work in parfor loops.
end


%
temp_A=[];  % data_1 container
care_A=[];
%
temp_B=[];  % data_2 container
care_B=[];
%
temp_C=[];  % data_3 container
care_C=[];



%% extract the points on the poincare section: by default the system
% status with zero velocities
% the period is around 6 seconds so keep the last ten cycles

%> -- filter the event points for zero velocity crossing
ind_zero_v                          = find(ieout==2);
len_ind_zero_v  = length(ind_zero_v);
if len_ind_zero_v
    trunc_num   = ceil(0.2*len_ind_zero_v);
    if len_ind_zero_v >1
        trunc_num    = min(trunc_num,8);
        trunct_index = len_ind_zero_v-trunc_num:len_ind_zero_v;
    else
        trunct_index = trunc_num;
    end
    care_A=[care_A,yeout(ind_zero_v(trunct_index),1)];
end

% avoid dismatch in dimension
logi=size(care_A);
if logi(1)==1;else;care_A=care_A';end

% adding section data after per calculation
temp_A=[temp_A,care_A];

%% y coordinate
% find the IC points near the fixed point
len_yeout0  = size(yeout0,1);
if len_yeout0
    trunc_num   = ceil(0.2*len_yeout0);
    if len_yeout0 >1
        trunct_index = len_yeout0-trunc_num:len_yeout0;
    else
        trunct_index = trunc_num;
    end
    care_B=[care_B;yeout0(trunc_num,2)];
    care_C=[care_C;yeout0(trunc_num,3)];
    
end
% avoid dismatch in dimension
logi=size(care_B);
if logi(1)==1;else;care_B=care_B'; end
% adding section data after per calculus
temp_B=[temp_B,care_B];

%% z coordinate
% avoid dismatch in dimension
logi=size(care_C);
if logi(1)==1; else; care_C=care_C'; end
% adding section data after per calculus
temp_C=[temp_C,care_C];

%% A collection
%
temp_A=[p_care,temp_A];
A_length_data(i)=length(temp_A);
A_diogram_matrix=[A_diogram_matrix,temp_A];
%

%% B collection
temp_B=[p_care,temp_B];
B_length_data(i)=length(temp_B);
B_diogram_matrix=[B_diogram_matrix,temp_B];

%% C collection
% get rid of duplicate and put velocity value as the first element
temp_C=[p_care,temp_C];
C_length_data(i)=length(temp_C);
C_diogram_matrix=[C_diogram_matrix,temp_C];


temp_run_info           = [p_care, run_info];
run_info_collection     = [run_info_collection,temp_run_info];

toc
figure; 
hold on
plot(yout(:,2), yout(:,1))
plot(yeout(:,2), yeout(:,1),'ro')
