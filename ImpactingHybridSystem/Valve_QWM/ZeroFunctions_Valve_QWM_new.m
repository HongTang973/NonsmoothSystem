function [prob_det] = ZeroFunctions_Valve_QWM (prob)

% par = {'a', 'r', 'kappa','L', 'T'}

%%  define the linear part & get the matrix A
% choice          = [1, 0, 0, 0, 1, 2,0];
% construct the default valve class without input or tune parameter by
% giving specific initialized p then pass on to the constructor
% par = [343.14, 0.78, 0.05 , 0.5,  T]';
% keys = {'a', 'rx', 'kappa','L', 'pd'};
keys = prob.keys;
par  = prob.par;
par_struct = struct_construction(keys, par);

p = nd_QWM_N1_gas(par_struct);
% clarify the equilibrium and the re-construnction order
y1 = p.x_max/p.x_ref;
%
p.y0 = [y1;    0;  p.delta + y1;    0;  0];
p.call_type = 'eq2jacobian';
% reconstruct the class according to the demand
if strcmp(p.call_type, 'eq2jacobian') % in constuction process, this loop
    % won't be activated;
    % in this case, the original structure of model is being passed on with
    % given state of the system
    P_  = (p.y0(3)+1)*p.p0;
    p.rho = p.get_rho(P_, p.T);
end
p = nd_QWM_N1_gas(p);
A = get_jacobian_at_steady_nd(p,p.y0);

% composed map Jacoi matrix and its eigenvalues
C  = [1 0 0 0 0];
% vx_= @(t,y)  C*A*(y-state0); 
B= -(1+p.rx)*[0;  1; 0 ; 0; 0];
% define initial condition: grid on the positive set
MW = B*C*A;
R = eye(length(C))+MW;
% sign_V = @(Y) sign(C*A_1*Y);
% %
A_s = (eye(5)-(B*C*A)/(C*A*B))*A;
pd = p.pd;
%
Mono_A = expm(A*pd);

M = R*Mono_A - eye(length(C));

%
prob_det = det(M);

end