clc
close all
%> initial guess for the SN point value: from the 
par     = [-0.1000 + 0.2000i, -0.1000 - 0.2000i, -0.5 , 1.781916719621720+0.000010259409,  1.6  , 5.847542373676547 ];
mu      = 0.00001;
eta     = -0.00001;
par_new = par;
par_new(3:4)    = par_new(3:4) + 0*[mu, eta];
[A,B,C,R]       = par2NForm(par_new);
%> initial rough approximation of the roots
[roots, fvals]  = get_roots_Det(A,R,C, [5 7],0);
abs(roots(2) - roots(1))
%> update the T parameter in the par list
par(6)          = roots(1);

%> check if the tangency point is the zero value point
% [roots, fvals]  = get_roots_Det(A,R,C, roots + [-0.01 0.01] , 1 );
% [fval, diff_det] = root_tangency(A,R,C, roots + [-0.01 0.01])
% [ diff_det_1] = root2tangency(A,R,C,roots(1))
% [ diff_det_2] = root2tangency(A,R,C,roots(2))

% 
%> find the more accurate value of critical b2 by getting the tangency
%> point to approach the zero value function
index               =  4;
%> this will give you the tangency point around the approximated roots
[fval_0, T_at_max_0]    =  tangency_point_in_range(roots + [-0.01 0.01], par,index, par(index))
[ diff_det_0]          = root2tangency(A,R,C,T_at_max_0)

%> use the dichotomy method to get the accurate critical b2 value
b2_span = [0.99*par(index) 1.01*par(index)];
b2det  = @(b2) tangency_point_in_range( roots + [-0.01 0.01], par, index, b2);
% b2det(b2_span(1))
% b2det(b2_span(2))
plot_flag = 0;
[b2_crit] = dichotomy(b2det, b2_span(1), b2_span(2), 1e-12, plot_flag);

% b2det(x)

[fval, T_at_max]          = tangency_point_in_range(roots + [-0.01 0.01], par,index, b2_crit);
[ tangency_at_tmax]       = root2tangency(A,R,C,T_at_max);

%> update the par0 with more accurate values
par_deli                  = par;
par_deli(index)           = b2_crit;
par_deli(end)             = T_at_max;
[A,B,C,R]                 = par2NForm(par_deli);
% [roots, fvals]  = get_roots_Det(A,R,C, roots + [-0.01 0.01] , 1 );
% par_deli = [-0.1 + 0.2i	-0.1 - 0.2i	-0.5 1.78192697901049 	1.6	5.84686867166138 ]
% n_tol = norm([fval, diff_det])

%> the zero function condition
% Det         = @(A,R,C,T) det(R*expm(A*T) - eye(length(C)));

%> to get the jacoian of the zero funciotnos at the SN point regarding the
%> parameter
Det(A,R,C,T_at_max )
root2tangency(A,R,C,  T_at_max )
indexes = [3, 4, 6];
[t_,Jacob]= get_jacobian(@F_zero, par_deli, indexes) %> & done
%> use the CO method to track the solution branch

%> pass the par_deli as the initial point for the continuation in the 2D
%> parameter plane
SN_codim2_plane;




%> 
function F = F_zero(par)

        [A,B,C,R]       = par2NForm(par);
        F               = [Det(A,R,C,par(6) );
            root2tangency(A,R,C,  par(6) )];

end

%

%> zero function condition
function f = Det(A,R,C,T)
f = det(R*expm(A*T) - eye(length(C)));

end

%> the tangency condition
function [ diff_det] = root2tangency(A,R,C,T)
Det         = @(t) det(R*expm(A*t) - eye(length(C)));
init_delta  = 1e-3;
min_delta   = 1e-12;
dt          = init_delta;
d_tmp           = 10;
iter            = 0;
while dt > min_delta
    d_new      = ( Det(T+dt) - Det(T-dt))/2/dt;
    tol        = abs(d_new - d_tmp);
    dt         = dt /2;
    d_tmp      = d_new;
    iter       = iter + 1;
    if tol < 1e-9
        break;
    end
    
end
diff_det = d_new;
end

%> derivative for single variable function
function varargout = tangency_point_in_range(tspan, par,index, var)
par(index)  = var;
[A,B,C,R]   = par2NForm(par);
fs          = 1e6;
delta       = 1/fs;
T           = tspan(1):delta:tspan(2);
%
Det         = @(t) det(R*expm(A*t) - eye(length(C)));

%> initialize the container of roots
% roots_ind   = [];
%> set the default sign as positive

fval_list   = zeros(1, length(T));
%> calculate the value of the det function and detect the sign change
for i=1:length(T)
    fval_list(i) = Det(T(i));
end

[fval, ind]  = max(fval_list);
%> output the results at the steady status with given valve lift
if nargout >= 0; varargout{1} =  fval;       end
if nargout >= 2; varargout{2} =  T(ind);            end
end
