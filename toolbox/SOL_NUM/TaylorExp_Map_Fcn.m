%> This is a function called to get the local Taylor expansion of a smooth
%> function/map to find the normal form of the bifurcation via projection
%> onto the center manifold
%> the Jcaobian A has been computed first
%> Fcn is given as a function of state variables and parameters
%>
function [B,C] = TaylorExp_Map_Fcn(prob)
T_p     = prob.sys_vec(prob.sys_vec_index.Tp_index);
sys_par = prob.sys_vec(prob.sys_vec_index.sys_par_index);
IC      = prob.sys_vec(prob.sys_vec_index.IC_index);
T       = prob.sys_vec(prob.sys_vec_index.T_simu_index);
%>
Fcn     = prob.Fcn_MaptobeExpanded;
x_p     = prob.Fcn_Map_fp;
%> calling method of this function is Fcn(IC, sys_par)
dim_sys  = prob.Fcn_Map_dim;
dim_par  = prob.Fcn_Map_dim_par;
%>
B = zeros(dim_sys,dim_sys^2);
C = zeros(dim_sys, dim_sys^3, dim_sys^3);
%>
%% > evaluate the norm form coefficients -- numerical approximation
init_delta  = norm(x_p)/10;
min_delta   = 1e-12;
delta       = init_delta;
for i = 1:dim_sys 
    for j = 1:dim_sys
        n_row = (i-1)*dim_sys + j;
        %> if i ~= j: (F(i+1,j+1) -F(i-1,j+1) -F(i+1,j-1)+F(i-1,j-1))/4/h^2
        %> if i  = j: (F(i+2) - 2*F(i) +F(i-2))/4/h^2
        i_vec = zeros(dim_sys, 1); i_vec(i) = 1;
        j_vec = zeros(dim_sys, 1); j_vec(i) = 1;
        B_ij_temp = ones(dim_sys,1);
        L_iter    = 0;
        while delta > min_delta
            %> approximate the second derivative f'' = [f(x + h) + f(x-h) - 2*f(x)]/h^2;
            B_ij_new = ( Fcn(x_p + (i_vec+j_vec)*delta) - Fcn(x_p + (j_vec - i_vec)*delta) ...
                - Fcn(x_p + (i_vec - j_vec)*delta) + Fcn(x_p - (j_vec + i_vec)*delta))/4/delta^2;

            %> approximate the third derivative  f'''(x) = [ f(x+2h) - f(x-2h) -2[f(x+h) - f(x-h)]]/2h^3;
            %     https://math.stackexchange.com/questions/1301769/...
            %     approximation-formula-for-third-derivative-is-my-approach-right

            

            %> check the convergence
            tol_1 = real(norm( B_ij_new - B_ij_temp ))/norm(B_ij_temp);
            
            %> update
            B_ij_temp  = B_ij_new;
            delta       = delta /2;
            %> update the iteration number
            L_iter = L_iter + 1;
            

            if tol_1 <= 1e-5
                B(:,n_row) = B_ij_temp;
                break;
            end
                

        end
        
    end
end
end
