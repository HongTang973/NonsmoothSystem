function [prob_det] = ZeroFunctions_3D_node (par)

% par = {'lambda_1', 'lambda_2', 'b2','b3', 'T'}
% par =  [-1, 0.2, 2.5, 0.5, 0.875926384431162];
%%  define the linear part & get the matrix A
lambda_3 = -1;
  
lambda_1 = par(1);

lambda_2 = par(2);


a1 = lambda_1 + lambda_2 + lambda_3;
a2 = -(lambda_1 * lambda_2 + lambda_2*lambda_3 + lambda_1*lambda_3);
a3 = lambda_1 * lambda_2 * lambda_3;
% 
b2 = par(3);
b3 =par(4);

A = [a1 1 0;a2 0 1; a3 0 0];
[V1,D1]=eig(A);
EA = @(T) real(V1*diag(exp(diag(D1)*T))*inv(V1));

B = [0; b2; b3];


C = [1,0,0];

R = eye(length(C)) - B*C*A;


% sign_V = @(Y) sign(C*A_1*Y);
% %
% A_s = (eye(8)-(B*C*A)/(C*A*B))*A;

T = par(5);

%
Mono_A = expm(A*T);
% Mono_A = EA(T);

M = R*Mono_A - eye(length(C));

%
prob_det = det(M);

end