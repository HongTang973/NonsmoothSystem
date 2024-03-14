%> test the composed_map function
SN_point_admis  = [-0.1 + 0.2i, -0.1 - 0.2i, -0.5 , 1.781916719621720+0.000010259409,  1.6  , 5.846837671661376];
equi_type       = 1;
%>
[A,B,C,R, T_2_det]  = par2NForm_Lienard(SN_point_admis);
x_00                = IC_generator(SN_point_admis(6),R,A,C,equi_type);
fprintf('Check the SN points! \n')
[Mono_p,Salt_p]     = IC2Floque_Multipliers(SN_point_admis(6),x_00,R,A,C)

%>
indexes     = [ 3, 4, 6];
[t_,Jacob]  = get_jacobian(@SIADS23_SN_zero_Fcns, SN_point_admis, indexes) %> & done

fprintf('The slope at (0,0) is %g \n', t_(2)/t_(1))

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
prob.evol_seq                       = 'map_phi';
% prob.evol_seq                       = 'phi_map';
prob.par_2_flow_operator            = @SIADS23_3D_ODEs;
prob.par_2_ResetMap                 = @SIADS23_3D_ResetMap;
prob.efunc                          = @SIADS23_3D_H_x;
prob.OB_C                           = [1, 0, 0];
%>
% Flow_ode      =  prob.par_2_flow_operator(par);
% ResetMap      =  prob.par_2_ResetMap(par);


%% > ---- start the specific problem and do the analysis
% x_p     = x_00;
x_p     = R\x_00;
% T_p     = 5.846837671661376;
% sys_par = [-0.1 + 0.2i, -0.1 - 0.2i, -0.5 , 1.781916719621720+0.000010259409,  1.6 ]';
T_p     = SN_point_admis(6);
sys_par = SN_point_admis(1:5)';
IC      = x_p;
T       = 20*T_p;

sys_vec = [x_p; T_p; sys_par; IC; T];

%% > evaluate the jacobian of the nonlinear map
prob.call_type = 'jacobian';
[J_pMP,eigD_pMP, eigV_pMP, dT_dy] = General_IHS_Numerical_PMP_Jacobian(prob, sys_vec);

%% Taylor expansion
prob.sys_vec = sys_vec;
prob.Fcn_Map_fp = x_p(C<1);
prob.Fcn_Map_dim = 2;
prob.Fcn_Map_dim_par = 5;
prob.Fcn_MaptobeExpanded = @(IC) PMP_SIADS23(prob,IC,sys_par);
[BB,CC] = TaylorExp_Map_Fcn(prob);
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
prob.fixed_point.v = v;
prob.fixed_point.w = w;
prob.fixed_point.eigV = eigV_pMP;
% prob.IPC_num_allowed = 1;
[poly_coeff_CtrErr, par_coeff_CtrErr,bilinear_coeff_CtrErr] = General_IHS_PMP_sys_par_coeff(prob, sys_vec)

[poly_coeff_2, par_coeff_2,bilinear_coeff_2] = IHS_PMP_3D_exapnsion_in_Eigenspace(prob, sys_vec)

%% > define the specific problem by setting the flow and the reset map
%  > and event function
a0 = par_coeff_CtrErr(3,1); % -0.1498      -0.1313
b0 = par_coeff_CtrErr(4,1); % -0.0472     -0.0414
a1 = bilinear_coeff_CtrErr(3,1); % phi_map -3.0289 map_phi 3.4711
b1 = bilinear_coeff_CtrErr(4,1); % 0.4224  0.0055 
c  = poly_coeff_CtrErr(1,2);
d  = poly_coeff_CtrErr(1,3);


-a0/b0
