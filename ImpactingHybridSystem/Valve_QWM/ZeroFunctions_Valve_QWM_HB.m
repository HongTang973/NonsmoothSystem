function [prob_det] = ZeroFunctions_Valve_QWM_HB (par)

% par = {'a', 'r', 'kappa','L', 'T'}

%%  define the linear part & get the matrix A
% choice          = [1, 0, 0, 0, 1, 2,0];
% construct the default valve class without input or tune parameter by
% giving specific initialized p then pass on to the constructor
% par = [343.14, 0.78, 0.05 , 0.5,  T]';
keys            = {'a', 'rx', 'kappa','L','RT','xg'};

par_struct      = struct_construction(keys, par);

p               = nd_QWM_N1_gas(par_struct);

    
% clarify the equilibrium and the re-construnction order
RT              = p.RT/100;
y1              = RT*p.x_max/p.x_ref;
%
p.y0            = [y1;    0;  p.delta + y1;    0;  0];
p.call_type     = 'eq2jacobian';
% reconstruct the class according to the demand
p               = nd_QWM_N1_gas(p);
% update the sonicvelocity if there is homotopy: now the gas density is 
%  updated due to its compressity
% due to potential gas & liquid mixture the net density and sonic velocity
% has to be modified according to the formula
rho_g   = p.rho;
rho_l   = 1000;
a_g     = p.a;
a_l     = 1300;
[rho_m, a_m, ~] =...
    get_rho_sonicvel_mixture(p.xg,rho_g,rho_l, a_g, a_l);
p.rho   = rho_m;
p.a     = a_m;
% update the dimensionless parameters
p       = nd_QWM_N1_gas(p);
A               = get_jacobian_at_steady_nd(p,p.y0);
[~,D]           = eig(A);
d = diag(D);
d(imag(d) == 0) =[];
prob_det        = prod(real(d));
end