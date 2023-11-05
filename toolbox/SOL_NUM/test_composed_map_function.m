%> test the composed_map function
clc
close all
%> so the input should be the problem and the parameter vector
PD_point_admis  = [0.1; 1; 0.3; 0.653363855944258; 0.788; 5.33165805420573];
[A,B,C,R, T_2_det] = prob.par2prob(PD_point_admis);
[~,~,~,x_00,~,~]    = LCO_Det_search(PD_point_admis(6),R,A,C,-1);
fprintf('Check the PD points! \n')
[Mono_p,Salt_p]     = IC2Floque_Multipliers(PD_point_admis(6),x_00,R,A,C)



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
prob.par_2_flow_operator            = @Shilnikov_chaos_ODEs;
prob.par_2_ResetMap                 = @Shilnikov_chaos_ResetMap;
prob.efunc                          = @H_x;
prob.OB_C                           = [1, 0, 0];
%>
% Flow_ode      =  prob.par_2_flow_operator(par);
% ResetMap      =  prob.par_2_ResetMap(par);


%% > ---- start the specific problem and do the analysis
x_p     = [-1;0.454520087056196;-0.535806513058751];
T_p     = 5.33165805420573;
sys_par = [0.1; 1; 0.3; 0.653363855944258; 0.788];
IC      = x_p;
T       = 2*T_p;

sys_vec = [x_p; T_p; sys_par; IC; T];

Maped_state = General_IHS_P1_Composed_Map(prob, sys_vec)
%> 

%> whenever you define the function to form the flow operator and the
%> imapct map, you already define the index of the parameters, the
%> dimension of the system



function  Flow_ode = Shilnikov_chaos_ODEs(par)
        %>
        [A,B,C,R,T_2_det] = par2NForm_Shilnikov(par);
         Flow_ode         = @(t,y) A*y;
        %>
end

function  ResetMap = Shilnikov_chaos_ResetMap(par)
        %>
        [A,B,C,R] = par2NForm_Shilnikov(par);
         ResetMap        = @(y0) R*y0;
end
%>
function [value,isterminal,direction] = H_x(t,y)
 value = real(y(1) + 1);                %  detect switching point
 isterminal = 1;                              %  stop the integration
 direction  = -1;
 end