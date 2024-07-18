%>
clc
%
prob.Fcn_MaptobeExpanded = @PolyFunc_test_1;
prob.Fcn_Map_fp          = [0;0;0];
%> calling method of this function is Fcn(IC, sys_par)
prob.Fcn_Map_dim  = 3;
prob.Fcn_Map_dim_par = 0;
[B,C] = TaylorExp_Map_Fcn(prob)