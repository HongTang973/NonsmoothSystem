% monitoring functions to track the change of cared feature of system

function [monitor] = Obeservation_Funs_Valve_QWM(par)

% par = {'U', 'r', 'damp_1','damp_2', 'damp_3', 'T'}

% output: cared values as a output vector
% floquet multipliers // amplitude of outcoming velocity // norm of the
% subspace vector // signV // A_s loci track

%%  define the linear part & get the matrix A
% choice          = [1, 0, 0, 0, 1, 2,0];
% construct the default valve class without input or tune parameter by
% giving specific initialized p then pass on to the constructor
% par = [343.14, 0.78, 0.05 , 0.5,  T]';
keys        = {'a', 'rx', 'kappa','L', 'pd'};

par_struct  = struct_construction(keys, par);

p           = nd_QWM_N1_gas(par_struct);
% clarify the equilibrium and the re-construnction order
y1          = p.x_max/p.x_ref;
%
p.y0        = [y1;    0;  p.delta + y1;    0;  0];
p.call_type = 'eq2jacobian';
% reconstruct the class according to the demand
p           = nd_QWM_N1_gas(p);
A           = get_jacobian_at_steady_nd(p,p.y0);

% composed map Jacoi matrix and its eigenvalues
C           = [1 0 0 0 0];
% vx_= @(t,y)  C*A*(y-state0); 
B           = -(1+p.rx)*[0;  1; 0 ; 0; 0];
% define initial condition: grid on the positive set
MW          = B*C*A;
R           = eye(length(C))+MW;
% sign_V = @(Y) sign(C*A_1*Y);
% %
A_s         = (eye(5)-(B*C*A)/(C*A*B))*A;
pd          = p.pd;
%
Mono_A      = expm(A*pd);
Mono_RA     = R*Mono_A;
M           = R*Mono_A - eye(length(C));

% construct a transform matrix
Trans       = eye(length(C));
Trans(C>0,:)= []; % get a (n-1) x n tranform matrix
%
ROW_M       = C*M; % to check the vadality of the remaining equation
COL_M       = Trans*M(:,C>0);

%
equi_type   =1;
M(C>0,:)    =[];
M(:,C>0)    =[];

% detect the singular point
% COND = rcond(M);
Ans         = -M\COL_M*equi_type;

sub_norm    = norm(Ans);

% candidate of the starting point of LCO
x_00        = equi_type*C'+ Trans'*Ans;

ve          = C*A*x_00;


% sing_v = sign_V(x_00);


F_ls        = @(y) A*y;
x_0         = inv(R)*x_00;
correction  = ((F_ls(x_00)-R*F_ls(x_0))*C)/(C*F_ls(x_0));
P           = R + correction;

PA          = P*Mono_A;

[~,D1]      = eig(Mono_RA);
Mono_p      = diag(D1);

[~,D2]      = eig(PA);

E2          = diag(D2);

[norm_salt, ind] =sort(abs(E2));

Salt_p      = E2(ind);

[~,D3]      = eig(A_s);

D_s         = diag(D3);

% 

monitor     = [Mono_p;norm_salt;Salt_p;ve;sub_norm;D_s];


end