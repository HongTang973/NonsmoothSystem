%% build-up of the coco in the non-smooth impacting hybrid system
% clear
clc
close all
%% -------------  get the coco's directory ----
% original=pwd;
% valve_project_startup;
% cd(original);

%% initialize the tool and customized settings
prob=coco_prob();
% prob = coco_set(prob, 'ode', 'vectorized', true); % vectorized or not the Jacobian
prob = coco_set(prob,'cont','h0',0.001); % 0.0001 or 0.01 different results
prob = coco_set(prob,'cont','h_min',10^-4); %0.0001
prob = coco_set(prob,'cont','h_max',1); %1
prob = coco_set(prob,'cont','h_fac_min',0.1); % Minimum step size adaptation factor ,0.5 100
prob = coco_set(prob,'cont','h_fac_max',2); % Maximum step size adaptation factor, 2 35 (can do strange stuff if too low..like not escaping)
prob = coco_set(prob,'cont','almax', 45);   % Critical angle between successive tangent vectors (try changing this one), 35
prob = coco_set(prob, 'cont', 'NullItMX', 1);
prob = coco_set(prob, 'coll', 'NTST', 30); %mesh 30, 100, 10000
prob = coco_set(prob, 'coll', 'NCOL', 5); %mesh 5,7
prob = coco_set(prob, 'coll', 'TOL', 10^-2); %mesh
% prob = coco_set(prob, 'cont', 'PtMX', 50); % maximum number of continuation steps
prob = coco_set(prob, 'cont', 'ItMX', 500); % maximum number of iteration - this is the one that impose the number of iterations
prob = coco_set(prob, 'cont', 'NPR', 1); %iteration
prob = coco_set(prob, 'cont', 'NAdapt', 1); % aadaptive changes are made to the orbit discretization after each successful step of continuation
prob = coco_set(prob,'corr','ItMX',50); %maximum number of iteration - for convergence (10 default)
prob = coco_set(prob,'corr','SubItMX',4); % damping in Newton criterion
prob = coco_set(prob,'corr','TOL',1.00E-009); %Tolerance of Newton correction
prob = coco_set(prob,'corr','ResTOL',1.00E-009); %Converge criterion of norm of the residium

% ------------- step 1: the zero problem ---
% Determine where your m-file's folder is.
% folder = fileparts(which(Obeservation_Funs_Airfoil.m));
% folder = 'F:\onedrive\OneDrive - University of Bristol\Codes stall\3DOF_AIRFOIL_MODEL\Bifurcation_Analysis\BEB\PO\';
% % Add that folder plus all subfolders to the path.
% addpath(genpath(folder));

%% continuation of the equlibrium branch
%  ----------   adding zero functions    ----------------------- %
% the command is as :
% ### coco_add_func(proc, fid, fhan, data, 'zero', 'uidx', uidx, 'u0', u0)
% ---fid: contains a function identifier
% ---fhan: contain a function handle
% ---data: contain the function data structure
% ---uidx: contain the index set
% ---u0: contain the vector u0
% ---'zero': function type

% ### coco_read_solution(fid, runid, labelID)

%> CO
% par                 = [-0.1000 + 0.2000i, -0.1000 - 0.2000i, -0.5 , 1.8,  1.6, 2.1314]';
% keys                = {'lambda_1', 'lambda_2', 'lambda_3','b2','b3','T'};
% par                 = [-0.1 , 0.2, -0.5 , 1.8,  1.6, 4.791524956512451];
par                 = [-0.1;0.0349756785977561;-0.5;1.8;1.6;7.51908282443018];
par = [-0.0895;
    0.0350;
   -0.5000;
    1.8000;
    1.6000;
   10.0000];
keys                = {'alpha', 'beta', 'lambda_3','b2','b3','T'};

%> pass the initial condition to the COCO
u0                  = par;

%> initialize the prob with customized settings

%% > embed the valve problem into the coco continuation problem

%> construct the problem from COCO built-in function
prob = coco_add_func(prob, 'Det_func', @COCO_3D_LCO_IHS, [], 'zero', 'u0', u0);

p_span   = [-5 5];
X_span   = [0 10];
% define continuation parameters
PNM             = keys;
PIdx1           = [1;2;3;4;5;6];
prob            = coco_add_pars(prob, '' ,PIdx1, PNM(PIdx1),'inactive');
bd0             = coco(prob, '3D_PD_case1', [], 1, {'b3','T'},  {p_span,X_span});

%> continue the case when chi is on

% prob1                           = prob;
% % prob1.COCO_Valve.prob           =  p1;
% bd1 = coco(prob1, '3D_case_b3_T', [], 1, {'b3','T'},  {p_span,X_span});
% % 
% %> continue the case when chi/Lambda  on
% prob2                           = prob;
% 
% bd2 = coco(prob2, '3D_case_lambda3_T', [], 1, {'lambda_3','T'},  {[-1 1],X_span});

figure; clf
haxes1 = axes;
thm = struct('special', {{'BP','FP'}});
coco_plot_bd(thm, '3D_PD_case1', 'b3', 'T')
hold on
% grid on
% ylim([0 0.8])
% xlim([0 120])
% figure; clf
% haxes2 = axes;
% thm = struct('special', {{'BP','FP'}});
% coco_plot_bd(thm, '3D_case_b3_T', 'b3', 'T')
% hold on
% 
% figure; clf
% haxes3 = axes;
% thm = struct('special', {{'BP','FP'}});
% coco_plot_bd(thm, '3D_case_lambda3_T', 'lambda_3', 'T')
% hold on
%>


%> read the fold bifurcation point
b2_list_bd0     = coco_bd_col(bd0, {'b2'});
T_list_bd0      = coco_bd_col(bd0, {'T'});
FPlab_bd0       = coco_bd_labs(bd0, 'FP');
lab_list_bd0            = coco_bd_col(bd0, 'LAB');
bull_ismember_bd0       = ismember(lab_list_bd0,FPlab_bd0);
FP_index_bd0            = find(bull_ismember_bd0);
b2_at_FP_bd0            = real(b2_list_bd0(FP_index_bd0));
T_at_FP_bd0             = real(T_list_bd0(FP_index_bd0));
%> so the first special codim-2 point on b2 X T plane
% keys                = {'lambda_1', 'lambda_2', 'lambda_3','b2','b3','T'};
% [-0.1000 + 0.2000i, -0.1000 - 0.2000i, -0.5 , 1.781916719621720,  1.6  , 5.847542373676547 ]';
%> 


%> read the fold bifurcation point
b3_list_bd1     = coco_bd_col(bd1, {'b3'});
T_list_bd1      = coco_bd_col(bd1, {'T'});
FPlab_bd1       = coco_bd_labs(bd1, 'FP');
lab_list_bd1            = coco_bd_col(bd1, 'LAB');
bull_ismember_bd1       = ismember(lab_list_bd1,FPlab_bd1);
FP_index_bd1            = find(bull_ismember_bd1);
b3_at_FP_bd1            = real(b3_list_bd1(FP_index_bd1));
T_at_FP_bd1             = real(T_list_bd1(FP_index_bd1));




%> so the first special codim-2 point on b3 X T plane
% keys                = {'lambda_1', 'lambda_2', 'lambda_3','b2','b3','T'};
% [-0.1000 + 0.2000i, -0.1000 - 0.2000i, -0.5 , 1.8,  1.591103794399196  , 5.698285027182902 ]';
%> 
%>

% lambda3_list_bd2         = coco_bd_col(bd2, {'lambda_3'});
% T_list_bd2               = coco_bd_col(bd2, {'T'});
% FPlab_bd2                = coco_bd_labs(bd2, 'FP');
% lab_list_bd2             = coco_bd_col(bd2, 'LAB');
% bull_ismember_bd2        = ismember(lab_list_bd2,FPlab_bd2);
% 
% FP_index_bd2             = find(bull_ismember_bd2 );
% lambda3_at_FP_bd2        = real(lambda3_list_bd2(FP_index_bd2));
% T_at_FP_bd2              = real(T_list_bd2(FP_index_bd2));

%> so the first special codim-2 point on lambda_3 X T plane
% keys                = {'lambda_1', 'lambda_2', 'lambda_3','b2','b3','T'};
% [-0.1000 + 0.2000i, -0.1000 - 0.2000i, -0.504103038128386 , 1.8,  1.6  , 5.676046546984852 ]';


