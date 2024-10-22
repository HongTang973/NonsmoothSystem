%> the code will record the 3D abstract model states on 3 chosen poincare
%> sections
%%  num for random initial conditions
num = 1; %< we give some initial condtions from pre-found P1 LCO
% define par loop parameter list
init_indx  = [];
par_1_list = [];
for i=1:num
    par_1_list=[par_1_list,p1_care_list];
end
% initcond
for j=1:num
    init_indx=[init_indx,j*ones(1,length(p1_care_list))];
end

%> initialize the buffer
%> initialize the variables
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
parfor i = 1:length(par_1_list)
    try
        p_care  = par_1_list(i);
        disp(p_care)
        %> choose one point with fixed eta
        point_1        = Par_ref;
        point_1(index) = point_1(index) + p_care*ds;

        %
        prob = full_BEB_SN_IHS_3D_par2prob_NL_David(point_1);
        %
        prob.v_th = nordmark_vth;
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
             care_A     = [care_A,yeout(ind_zero_v,1)];
            % trunc_num   = ceil(0.2*len_ind_zero_v);
            % if len_ind_zero_v >1
            %     trunc_num    = min(trunc_num,8);
            %     trunct_index = len_ind_zero_v-trunc_num:len_ind_zero_v;
            % else
            %     trunct_index = trunc_num;
            % end
            % care_A=[care_A,yeout(ind_zero_v(trunct_index),1)];
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
            care_B = [care_B;yeout0(:,2)];
            care_C = [care_C;yeout0(:,3)];
            % trunc_num   = ceil(0.2*len_yeout0);
            % if len_yeout0 >1
            %     trunct_index = len_yeout0-trunc_num:len_yeout0;
            % else
            %     trunct_index = trunc_num;
            % end
            % care_B=[care_B;yeout0(trunc_num,2)];
            % care_C=[care_C;yeout0(trunc_num,3)];

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

    catch ME
        disp(['error at',num2str(p_care)])
        rethrow(ME)
    end
end
toc

