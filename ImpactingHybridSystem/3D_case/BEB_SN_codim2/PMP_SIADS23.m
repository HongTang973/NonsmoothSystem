function f = PMP_SIADS23(prob,IC,sys_par)
x_p     = prob.sys_vec(prob.sys_vec_index.xp_index);
T_p     = prob.sys_vec(prob.sys_vec_index.Tp_index);
T       = prob.sys_vec(prob.sys_vec_index.T_simu_index);
C       = prob.OB_C;
%
ind         = C<1;
IC_flow     = x_p;
IC_flow(ind)= IC; % n-1 dimensional vector
%
sys_vec = [x_p; T_p; sys_par; IC_flow; T];

f = General_IHS_P1_Composed_Map(prob, sys_vec);

f(end) = [];
f(C>0) = [];
end 