function [poly_coeff_CtrErr, par_coeff_CtrErr,bilinear_coeff_CtrErr] = General_IHS_PMP_sys_par_coeff(prob, sys_vec)
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
prob.call_type = 'Coeff_expansion';
%% > ------------------ OUTPUT ----------------------------------------- %
% > the coefficients of x^2  x^3
% > the coefficients of the par_1 par_2 ... par_n


% sys_vec = [x_p; T_p; sys_par; IC; T];
x_p     = sys_vec(prob.sys_vec_index.xp_index);
T_p     = sys_vec(prob.sys_vec_index.Tp_index);
sys_par = sys_vec(prob.sys_vec_index.sys_par_index);
IC      = sys_vec(prob.sys_vec_index.IC_index);
T       = sys_vec(prob.sys_vec_index.T_simu_index);
C       = prob.OB_C;

% dim_sys : the dimension of the system
dim_sys        = length(x_p);
dim_par        = length(sys_par);
N_loop         = 100;

v = prob.fixed_point.v;
w = prob.fixed_point.w;

%% > evaluate the norm form coefficients -- numerical approximation
init_delta  = norm(x_p)/10;
min_delta   = 1e-12;
ind         = find(C<1);
S           = zeros(dim_sys, dim_sys-1);
S(ind,:)    = eye(dim_sys-1);

%>  ---------------------------------------------------------------- %
L_iter      = 0;
lambda_tmp  = 0;
c_tmp       = 0;
d_tmp       = 0;
%>
tol_lambda  = 1;
tol_c       = 1;
tol_d       = 1;
delta       = init_delta;
%>
while delta > min_delta
    %  ------- Phi(x_p + delta) -------
    phi_p1 = General_IHS_P1_Composed_Map(prob, ...
        [x_p; T_p; sys_par; x_p + S*delta* v; T]);
    %  ------- Phi(x_p - delta) -------
    phi_m1 = General_IHS_P1_Composed_Map(prob, ...
        [x_p; T_p; sys_par; x_p - S*delta* v; T]);
    %  ------- Phi(x_p + 2*delta) -------
    phi_p2 = General_IHS_P1_Composed_Map(prob, ...
        [x_p; T_p; sys_par; x_p + 2*S*delta* v; T]);
    %  ------- Phi(x_p - 2*delta) -------
    phi_m2 = General_IHS_P1_Composed_Map(prob, ...
        [x_p; T_p; sys_par; x_p - 2*S*delta* v; T]);
    %> approximate the first derivative f' = [f(x + h) - f(x-h)]/2/h;
    lambda_new = w'*( phi_p1(ind) - phi_m1(ind))/2/delta;

    %> approximate the first derivative f'' = [f(x + h) + f(x-h) - 2*f(x)]/h^2;
    c_new      = w'*( phi_p1(ind) +  phi_m1(ind) - 2*x_p(ind))/delta^2;

    %> approximate the first derivative  f'''(x) = [ f(x+2h) - f(x-2h) -2[f(x+h) - f(x-h)]]/2h^3;
    %     https://math.stackexchange.com/questions/1301769/...
    %     approximation-formula-for-third-derivative-is-my-approach-right

    d_new      = w'*( phi_p2(ind) - phi_m2(ind) - 2*phi_p1(ind) + 2*phi_m1(ind) )/2/delta^3;

    %> check the convergence
    tol_1 = real(norm( lambda_new - lambda_tmp )/norm(lambda_new));
    tol_2 = real(norm( c_new - c_tmp )/norm(c_new));
    tol_3 = real(norm( d_new - d_tmp )/norm(d_new));
    %> update
    lambda_tmp  = lambda_new;
    c_tmp       = c_new;
    d_tmp       = d_new;
    delta       = delta /2;
    %> update the iteration number
    L_iter = L_iter + 1;
    %> check the tolerence
    if tol_3 < tol_d
        d     =  d_tmp;
        tol_d = tol_3;
    end
    %> check the tolerence
    if tol_2 < tol_c
        c     =  c_tmp;
        tol_c = tol_2;
    end
   
    %> check the tolerence
    if tol_1 < tol_lambda
        lambda     =  lambda_tmp;
        tol_lambda = tol_1;
    end
    
    if tol_lambda <= 1e-6 && tol_c <= 1e-6 && tol_d <= 1e-6
        break;
    end


end

% lambda = lambda_new;
%  c      = c_new;
%  d      = d_new;
poly_coeff_CtrErr      = [lambda, c/2, d/6; real(tol_lambda), real(tol_c), real(tol_d)];
fprintf('lambda converges to %g with %d times iteration and RelTol %g \n', lambda, L_iter, tol_lambda)


%% > evaluate the first derivative regarding the parameters of the system
init_delta  = 1e-2;
min_delta   = 1e-12;
par_coeff_CtrErr  = zeros(dim_par,6);

for j = 1:dim_par

    delta               =   init_delta;
    par_selection       = zeros(size(sys_par));
    par_selection(j)    = 1;
    FO_tmp              = 0;
    SO_tmp              = 0;
    Te_tmp              = 0;
    tol_FO              = 10;
    tol_SO              = 10;
    tol_Te              = 10;
    while delta > min_delta
        %  ------- Phi(sys_par + delta) -------
        phi_p1 = General_IHS_P1_Composed_Map(prob, ...
            [x_p; T_p; sys_par + delta*par_selection; x_p ; T]);
        %  ------- Phi(sys_par - delta) -------
        phi_m1 = General_IHS_P1_Composed_Map(prob, ...
            [x_p; T_p; sys_par - delta*par_selection; x_p ; T]);
        %      %  ------- Phi(sys_par + 2*delta) -------
        %     phi_p2 = General_IHS_P1_Composed_Map(prob, ...
        %         [x_p; T_p; sys_par; x_p + 2*S*delta* v; T]);
        %      %  ------- Phi(sys_par - 2*delta) -------
        %     phi_m2 = General_IHS_P1_Composed_Map(prob, ...
        %         [x_p; T_p; sys_par; x_p - 2*S*delta* v; T]);
        %> approximate the first derivative f' = [f(x + h) - f(x-h)]/2/h;
        FO_new = w'*( phi_p1(ind) - phi_m1(ind))/2/delta;

        %> approximate the first derivative f'' = [f(x + h) + f(x-h) - 2*f(x)]/h^2;
        SO_new      = w'*( phi_p1(ind) +  phi_m1(ind) - 2*x_p(ind))/delta^2;

        %> approximate the first derivative  f'''(x) = [ f(x+2h) - f(x-2h) -2[f(x+h) - f(x-h)]]/2h^3;
        %     https://math.stackexchange.com/questions/1301769/...
        %     approximation-formula-for-third-derivative-is-my-approach-right

        % d_new      = w'*( phi_p2(ind) - phi_m2(ind) - 2*phi_p1(ind) + 2*phi_m1(ind) )/2/delta^3;
        Te_new = ( phi_p1(end) - phi_m1(end))/2/delta;
        %> check the convergence
        if norm( FO_new) > 1e-3
            tol_1 = real(norm( FO_new - FO_tmp )/norm(FO_new));
        else
            tol_1 = real(norm( FO_new - FO_tmp ));
        end
        %
        if norm( SO_new) > 1e-3
            tol_2 = real(norm( SO_new - SO_tmp )/norm(SO_new));
        else
            tol_2 = real(norm( SO_new - SO_tmp ));
        end
        % 
        if norm( Te_new) > 1e-3
            tol_3 = real(norm( Te_new - Te_tmp )/norm(Te_new));
        else
            tol_3 = real(norm( Te_new - Te_tmp ));
        end
        %> update
        FO_tmp       = FO_new;
        SO_tmp       = SO_new;
        Te_tmp       = Te_new;
        % d_tmp       = d_new;
        delta       = delta /2;
        %> update the iteration number
        L_iter = L_iter + 1;
        
        %> check the tolerence
        if tol_3 <= tol_Te
            Te     =  Te_tmp;
            tol_Te = tol_3;
        end
        
        %> check the tolerence
        if tol_2 <= tol_SO
            SO     =  SO_tmp;
            tol_SO = tol_2;
        end
        %> check the tolerence
        if tol_1 <= tol_FO
            FO    = FO_tmp;
            tol_FO = tol_1;
        end

        %>
        if tol_FO <=1e-6 && tol_SO <=1e-6 && tol_Te <=1e-6
            break;
        end


    end
    par_coeff_CtrErr(j,1) = FO;
    par_coeff_CtrErr(j,2) = SO/2;
    par_coeff_CtrErr(j,3) = Te;
    par_coeff_CtrErr(j,4) = real(tol_FO);
    par_coeff_CtrErr(j,5) = real(tol_SO);
    par_coeff_CtrErr(j,6) = real(tol_Te);
end
    %% > get the bilinear derivative which should be a matrix : J*u*par_i
    init_delta  = 1e-2;
    min_delta   = 1e-12;
    bilinear_coeff_CtrErr  = zeros(dim_par,2);

    for j = 1:dim_par

        delta               =   init_delta;
        par_selection       = zeros(size(sys_par));
        par_selection(j)    = 1;
        FO_tmp              = 0;
        SO_tmp              = 0;
        tol_FO              = 10;
        tol_SO              = 10;
        while delta > min_delta
            %  ------- Phi(sys_par + delta, x_p + delta) -------
            phi_p1_p1 = General_IHS_P1_Composed_Map(prob, ...
                [x_p; T_p; sys_par + delta*par_selection; x_p + S*delta* v ; T]);
            %  ------- Phi(sys_par - delta, x_p - delta) -------
            phi_m1_m1 = General_IHS_P1_Composed_Map(prob, ...
                [x_p; T_p; sys_par - delta*par_selection; x_p - S*delta* v ; T]);
            %  ------- Phi(sys_par + delta, x_p - delta) -------
            phi_p1_m1 = General_IHS_P1_Composed_Map(prob, ...
                [x_p; T_p; sys_par + delta*par_selection; x_p - S*delta* v ; T]);
            %  ------- Phi(sys_par - delta, x_p + delta) -------
            phi_m1_p1 = General_IHS_P1_Composed_Map(prob, ...
                [x_p; T_p; sys_par - delta*par_selection; x_p + S*delta* v ; T]);
            %> approximate the first derivative f' = [f(x + h) - f(x-h)]/2/h;
            tmp         = ( phi_p1_p1(ind) +  phi_m1_m1(ind) - phi_p1_m1(ind) - phi_m1_p1(ind))/4/delta^2;
            FO_new      = w'*tmp;

            %> check the convergence
            if norm( FO_new) > 1e-3
                tol_1 = real(norm( FO_new - FO_tmp )/norm(FO_new));
            else
                tol_1 = real(norm( FO_new - FO_tmp ));
            end
         
            %> update
            FO_tmp       = FO_new;
            delta       = delta /2;
            %> update the iteration number
            L_iter = L_iter + 1;

         
            %> check the tolerence
            if tol_1 <= tol_FO
                FO    = FO_tmp;
                tol_FO = tol_1;
            end

            %>
            if tol_FO <=1e-6 
                break;
            end


        end
        bilinear_coeff_CtrErr(j,1) = FO;
        bilinear_coeff_CtrErr(j,2) = real(tol_FO);
    end

