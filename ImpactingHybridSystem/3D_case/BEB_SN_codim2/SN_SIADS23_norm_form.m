%> test the composed_map function
% clc
close all
% clear
%> so the input should be the problem and the parameter vector
SN_point_admis  = [-0.1 - 0.2i,-0.1 + 0.2i,-0.5,1.78192697900802,1.6,5.84691856641790 ];
equi_type       = 1;
%>
[A,B,C,R, T_2_det]  = par2NForm_Lienard(SN_point_admis);
[~,~,~,x_00,~,~]    = LCO_Det_search(SN_point_admis(6),R,A,C,equi_type);
fprintf('Check the PD points! \n')
[Mono_p,Salt_p]     = IC2Floque_Multipliers(SN_point_admis(6),x_00,R,A,C)



%% -> define the prob
dim_sys  = 3;
dim_par  = 5;
%>
xp_index =          1:1:dim_sys;
Tp_index =          xp_index(end) + 1;
sys_par_index =     [1:1:dim_par] + Tp_index(end);
IC_index      =     [1:1:dim_sys] + sys_par_index(end);
T_simu_index  = IC_index(end) +1;
%>  ----- define the index of the state variables and the system parameters
prob.sys_vec_index.xp_index         = xp_index;
prob.sys_vec_index.Tp_index         = Tp_index;
prob.sys_vec_index.sys_par_index    = sys_par_index ;
prob.sys_vec_index.IC_index         = IC_index;
prob.sys_vec_index.T_simu_index     = T_simu_index;

%> ----- define the functions to form the question
prob.par_2_flow_operator            = @SIADS23_3D_ODEs;
prob.par_2_ResetMap                 = @SIADS23_3D_ResetMap;
prob.efunc                          = @SIADS23_3D_H_x;
prob.OB_C                           = [1, 0, 0];
%>
% Flow_ode      =  prob.par_2_flow_operator(par);
% ResetMap      =  prob.par_2_ResetMap(par);


%% > ---- start the specific problem and do the analysis
x_p     = x_00;
T_p     = 5.846837671661376;
sys_par = [-0.1 + 0.2i, -0.1 - 0.2i, -0.5 , 1.781916719621720+0.000010259409,  1.6 ]';
IC      = x_p;
T       = 20*T_p;

sys_vec = [x_p; T_p; sys_par; IC; T];

%% > evaluate the jacobian of the nonlinear map
prob.call_type = 'jacobian';
[J_pMP,eigD_pMP, eigV_pMP, dT_dy] = General_IHS_Numerical_PMP_Jacobian(prob, sys_vec);

%
%> find the right & left eigenvector and normalize the left eigenvector
[eigLV_pMP, eigLD_pMP] = eig(J_pMP');
%> locate the index of the eigen vector corresponding to the -1 multiplier
D_eig = diag(eigD_pMP);
PD_ind = find(abs(D_eig -1) < 1e-3);
LD_eig = diag(eigLD_pMP);
PD_ind_L = find(abs(LD_eig -1) < 1e-3);
%> pick out the w/v
v = eigV_pMP(:,PD_ind);
w = eigLV_pMP(:,PD_ind_L);
w = w/norm(w'*v);

%> 
prob.PD_fixed_point.v = v;
prob.PD_fixed_point.w = w;
[poly_coeff_CtrErr, par_coeff_CtrErr] = General_IHS_PMP_sys_par_coeff(prob, sys_vec)



%% > define the specific problem by setting the flow and the reset map
%  > and event function