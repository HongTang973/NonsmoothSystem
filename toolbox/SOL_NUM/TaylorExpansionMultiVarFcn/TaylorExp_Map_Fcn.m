%> This is a function called to get the local Taylor expansion of a smooth
%> function/map to find the normal form of the bifurcation via projection
%> onto the center manifold
%> the Jcaobian A has been computed first
%> Fcn is given as a function of state variables and parameters
%>
function [B,C,B_err,C_err] = TaylorExp_Map_Fcn(prob)
%>
Fcn     = prob.Fcn_MaptobeExpanded;
x_p     = prob.Fcn_Map_fp;
%> calling method of this function is Fcn(IC, sys_par)
dim_sys  = prob.Fcn_Map_dim;
dim_par  = prob.Fcn_Map_dim_par;
%>
B = zeros(dim_sys,dim_sys^2);
C = zeros(dim_sys,dim_sys^3);
B_err = zeros(1,dim_sys^2);
C_err = zeros(1,dim_sys^3);
%>
%% > evaluate the norm form coefficients -- numerical approximation
init_delta  = min(1, 10^norm(x_p)/100);
min_delta   = 1e-12;

for i = 1:dim_sys
    for j = 1:dim_sys
        B_n_row = (i-1)*dim_sys + j;
        %> if i ~= j: (F(i+1,j+1) -F(i-1,j+1) -F(i+1,j-1)+F(i-1,j-1))/4/h^2
        %> if i  = j: (F(i+2) - 2*F(i) +F(i-2))/4/h^2
        i_vec = zeros(dim_sys, 1); i_vec(i) = 1;
        j_vec = zeros(dim_sys, 1); j_vec(j) = 1;
        B_ij_temp = ones(dim_sys,1);
        B_iter    = 0;
        tol_B     = 1;
        delta     = init_delta;
        while delta > min_delta
            %> approximate the second derivative f'' = [f(x + h) + f(x-h) - 2*f(x)]/h^2;
            B_ij_new = (  Fcn(x_p + ( i_vec + j_vec)*delta) ...
                - Fcn(x_p + (-i_vec + j_vec)*delta) ...
                - Fcn(x_p + ( i_vec - j_vec)*delta) ...
                + Fcn(x_p + (-j_vec - i_vec)*delta))...
                /4/delta^2;
            %> check the convergence
            tol_1 = real(norm( B_ij_new - B_ij_temp ))/norm(B_ij_temp);
            %> update
            B_ij_temp   = B_ij_new;
            delta       = delta /2;
            %> update the iteration number
            B_iter      = B_iter + 1;
            %> check the tolerence
            if tol_1 < tol_B
                B(:,B_n_row)   = B_ij_temp;
                tol_B          = tol_1;
                B_err(B_n_row) = tol_B;
            end
            
            if tol_1 <= 1e-6; break; end
            %
        end
        %
        for k = 1 : dim_sys
            %> approximate the third derivative  f'''(x) = [ f(x+2h) - f(x-2h) -2[f(x+h) - f(x-h)]]/2h^3;
            %     https://math.stackexchange.com/questions/1301769/...
            %     approximation-formula-for-third-derivative-is-my-approach-right

            C_n_row = (i-1)*dim_sys + j + (k-1)*dim_sys^2;
            %> if i ~= j: (F(i+1,j+1) -F(i-1,j+1) -F(i+1,j-1)+F(i-1,j-1))/4/h^2
            %> if i  = j: (F(i+2) - 2*F(i) +F(i-2))/4/h^2
            k_vec     = zeros(dim_sys, 1); k_vec(k) = 1;
            C_ijk_temp = ones(dim_sys,1);
            C_iter    = 0;
            tol_C     = 1;
            delta     = init_delta;
            while delta > min_delta
                %> approximate the second derivative f'' = [f(x + h) + f(x-h) - 2*f(x)]/h^2;
                C_ijk_new = ( Fcn(x_p + ( i_vec + j_vec + k_vec)*delta) ...
                    - Fcn(x_p + (-i_vec + j_vec + k_vec)*delta) ...
                    - Fcn(x_p + ( i_vec - j_vec + k_vec)*delta) ...
                    + Fcn(x_p + (-i_vec - j_vec + k_vec)*delta) ...
                    - Fcn(x_p + ( i_vec + j_vec - k_vec)*delta) ...
                    + Fcn(x_p + (-i_vec + j_vec - k_vec)*delta) ...
                    + Fcn(x_p + ( i_vec - j_vec - k_vec)*delta) ...
                    - Fcn(x_p + (-i_vec - j_vec - k_vec)*delta) ...
                    )/8/delta^3;
                %> check the convergence
                tol_2 = real(norm( C_ijk_new - C_ijk_temp ))/norm(C_ijk_temp);

                %> update
                C_ijk_temp  = C_ijk_new;
                delta       = delta /2;
                %> update the iteration number
                C_iter      = C_iter + 1;

                %> check the tolerence
                if tol_2 < tol_C
                    C(:,C_n_row)    = C_ijk_temp;
                    tol_C           = tol_2;
                    C_err(C_n_row)  = tol_C;
                end
                if tol_2 <= 1e-6; break; end
            end %> end the loop to find the third derivative
        end %> third index
    end %> second index
end %> first index
end
