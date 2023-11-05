function  [A,R,C]= Matrices_3D_impact(lambda_1,lambda_2,lambda_3,b2,b3,C)

% lambda_3 = -1;
%   
% lambda_1 = -1;
% 
% lambda_2 = 0.2;


a1 = lambda_1 + lambda_2 + lambda_3;
a2 = -(lambda_1 * lambda_2 + lambda_2*lambda_3 + lambda_1*lambda_3);
a3 = lambda_1 * lambda_2 * lambda_3;
% 
% b2 = 2.5; b3 =1;

A = [a1 1 0;a2 0 1; a3 0 0];

B = [0; b2; b3];


% C = [1,0,0];

R = eye(length(C)) - B*C*A;