function [J_pMP,eigD_pMP, eigV_pMP, dT_dy] = General_IHS_Numerical_PMP_Jacobian(prob, sys_vec)
% For a general Hybrid system with a flow region and the reset map on the
% discontinuity surface, if there is a LCO with
% x_0: incoming point
% x_00: outcoming point
% x_p : point on the poincare section
% T: the period of this orbit
% F_ls: function handle of the flow
% RMap: the reset map when the flow hit the boundary
% EventFun: the function to detect the event of hitting the boundary
% next two variables to plot out phase, last integers to specify the
% figure number to open
% *********************************************************%
%---------------------------AIM----------------------------%
% Calculate the Poincare map numerically around the LCO to check %
% the stability of the LCO                                 %
% *********************************************************%

% sys_vec = [x_p; T_p; sys_par; IC; T];
x_p     = sys_vec(prob.sys_vec_index.xp_index);
T_p     = sys_vec(prob.sys_vec_index.Tp_index);
sys_par = sys_vec(prob.sys_vec_index.sys_par_index);
IC      = sys_vec(prob.sys_vec_index.IC_index);
T       = sys_vec(prob.sys_vec_index.T_simu_index);
C       = prob.C;

% n_s : the dimension of the system
n_s         = length(x_p);
% define the num of loop and store the LPE
N_loop      = 100;
init_delta  = 1; 
min_delta   = 1e-12;

%> the tmp of Jacobian
ind         = find(C<1);
J_pMP       = zeros(n_s-1,n_s-1);
dT_dy       = zeros(1,n_s-1);
%> 
Jac_temp    = J_pMP;

%>
for i = 1:n_s-1 
    delta = init_delta;
    %>
    vec_selection = zeros(n_s,1);
    vec_selection(ind(i)) = 1;
    while delta > min_delta
        phi_p1 = General_IHS_P1_Composed_Map(prob, ...
            [x_p; T_p; sys_par; x_p + delta* vec_selection; T]);
        phi_m1 = General_IHS_P1_Composed_Map(prob, ...
            [x_p; T_p; sys_par; x_p - delta* vec_selection; T]);
        %> 
        Jac_new     = ( phi_p1(ind)- phi_m1(ind))/2/delta;
        dT_dy_new   = ( phi_p1(end)- phi_m1(end) )   /2/delta;
        
        %> check if the Jac converges
        tol = norm( Jac_new - Jac_temp(:,i) )/norm(Jac_new);
        %> update
        Jac_temp(:,i) = Jac_new;
        dT_dy (i)     = dT_dy_new;
        delta = delta / 2;
        %>
        if tol < 1e-6
            break;
        end
        
    end
    
end
% 
 J_pMP = Jac_temp;
%
[eigV_pMP, eigD_pMP] = eig(J_pMP);