%> ---------------- BRUTE FORCE DIAGRAM SIMULATION ----------------- %
clear
clc
close all

%% > initialization for the parallel computing
Initialization_for_Par;

%> include the path needed for running
run('../startup.m')
sprintf('The corresponding paths have been included!')
%> initialize the structure of the gas object
p = get_initialized_p;
%> define the initial condition as a perturbation around the equilibrium
p.Qm = 1.0;
% F_quadratic_nd_QWM(0,p.y0,p)
% IC         = p.y0 + [0 -0.2 0.1 0 0]';
% F_quadratic_nd_QWM(0, IC, p)
tspan      = [0 100];
fs         = 40;
p.gap      = p.x_max/p.x_ref;


%%  define the parameter range interested --- 1D  parameter sweeping

p1_up_limit  = 1.5;
p1_bot_limit = 0.5;
delta_1   = 0.01;
p1_care_list=p1_bot_limit:delta_1:p1_up_limit;

% num for random initial conditions
num = 7;
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

%> initialize the 
beta_diogram_matrix=[]; % vector to keep result
beta_length_data=[];    % record length for per result in order to recover
% data from vector to matrix
%
B_diogram_matrix=[];
B_length_data=[];
%
C_diogram_matrix=[];
C_length_data=[];

%> start the computation using paralleledcomputation
tic
parfor i = 1:length(par_1_list)
    try
        L_care=par_1_list(i);
        disp(L_care)


        %> ----- reform the problem using the defined parameter  ---------%

        p_temp   = p;
        p_temp.L = L_care;
        p_temp   = nd_QWM_N1_gas(p_temp)

            
        %> specify the int function handle other than the default ones
        p_temp.func_han = @F_LN_nd_QWM;
        %
        temp_A=[];  % data_1 container
        care_A=[];
        %
        temp_B=[];  % data_2 container
        care_B=[];
        %
        temp_C=[];  % data_3 container
        care_C=[];

        n_init=init_indx(i);

        
        IC         = p_temp.y0 + [0 -0.2*(n_init-1)/(num-1) 0.01 0 0]';
        %> use the full model to run the simulation
        [tout,yout,~,~,~,~,~]=...
            Valve_int_ode_QWM(p_temp,IC,tspan,fs);

        %% extract the points on the poincare section: by default the system 
        % can is integratable and the zero velocity points can be captured
        % by the built-in event detecting function

        if tout(end)>tspan(2)-1
            % successfully integrate to the end
            idx1=find(tout>0.8*tout(end));
        end

        %> get the local minimum and maximum
        data=yout(idx1,1)';
        IndMin=find(diff(sign(diff(data)))>0)+1;   %local minumum index
        IndMax=find(diff(sign(diff(data)))<0)+1;   %local maximum index
        care_A=[care_A,data(IndMin)];
        care_A=[care_A,data(IndMax)];

        if isempty(IndMin) && isempty(IndMax)
            % means that get into equilibrium
            care_A=[care_A,data(end)];
        end
        % set tol for abrevity
        care_A(abs(care_A)<1e-12)=0;
        % delete same  elements
        care_A=unique(care_A);
        % avoid dismatch in dimension
        logi=size(care_A);
        if logi(1)==1
        else
            care_A=care_A';
        end
        % adding section data after per calculus
        temp_A=[temp_A,care_A];

        %% pitch motion
        % find v==0 points with set tol ----poincare section method I
        C_idx2=abs(yout(:,5))<1e-6;
        care_C=yout(C_idx2,2);
        % poincare section method II
        data=yout(:,2);
        IndMin=find(diff(sign(diff(data)))>0)+1;   %获得局部最小值的位置
        IndMax=find(diff(sign(diff(data)))<0)+1;   %获得局部最大值的位置
        care_C=[care_C;data(IndMin)];
        care_C=[care_C;data(IndMax)];
        % set tol for abrevity
        care_C(abs(care_C)<1e-12)=0;
        % delete same zero elements
        care_C=unique(care_C);
        % avoid dismatch in dimension
        logi=size(care_C);
        if logi(1)==1
        else
            care_C=care_C';
        end
        % adding section data after per calculus
        temp_C=[temp_C,care_C];

        %% B motion
        % find v==0 points with set tol ----poincare section method I
        B_idx2=abs(yout(:,4))<1e-6;
        care_B=yout(B_idx2,1);
        % poincare section method II
        data=yout(:,1);
        IndMin=find(diff(sign(diff(data)))>0)+1;   %获得局部最小值的位置
        IndMax=find(diff(sign(diff(data)))<0)+1;   %获得局部最大值的位置
        care_B=[care_B;data(IndMin)];
        care_B=[care_B;data(IndMax)];
        % set tol for abrevity
        care_B(abs(care_B)<1e-12)=0;
        % delete same zero elements
        care_B=unique(care_B);
        % avoid dismatch in dimension
        logi=size(care_B);
        if logi(1)==1
        else
            care_B=care_B';
        end
        % adding section data after per calculus
        temp_B=[temp_B,care_B];
        
        %% beta collection
        % get rid of duplicate and put velocity value as the first element
        
        temp_A=sort(temp_A,'descend');
        front_serial=temp_A(1:end-1);
        latter_serial=temp_A(2:end);
        adjacent_dis=abs(front_serial-latter_serial);
        temp_A=[temp_A(1),latter_serial(adjacent_dis>1e-6)];
        
        %
        temp_A=[L_care,temp_A];
        beta_length_data(i)=length(temp_A);
        beta_diogram_matrix=[beta_diogram_matrix,temp_A];
        %
        
        %% C collection
        % get rid of duplicate and put velocity value as the first element
        temp_C=sort(temp_C,'descend');
        front_serial=temp_C(1:end-1);
        latter_serial=temp_C(2:end);
        adjacent_dis=abs(front_serial-latter_serial);
        temp_C=[temp_C(1),latter_serial(adjacent_dis>1e-6)];
        %
        temp_C=[L_care,temp_C];
        C_length_data(i)=length(temp_C);
        C_diogram_matrix=[C_diogram_matrix,temp_C];
        
        %% B collection
        % get rid of duplicate and put velocity value as the first element
        temp_B=sort(temp_B,'descend');
        front_serial=temp_B(1:end-1);
        latter_serial=temp_B(2:end);
        adjacent_dis=abs(front_serial-latter_serial);
        temp_B=[temp_B(1),latter_serial(adjacent_dis>1e-6)];
        %
        temp_B=[L_care,temp_B];
        B_length_data(i)=length(temp_B);
        B_diogram_matrix=[B_diogram_matrix,temp_B];
        
    catch ME
        disp(['error at',num2str(L_care)])
    end
end
toc
save(sprintf('L_%s_.mat',date))
