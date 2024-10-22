function  [A,B,C,R,A_s]= Matrices_3D_IHS_Gform(p,flag)

% built-in default poincare section place
phi = p.phi;
theta = p.theta;
e1 = [cos(theta)*sin(phi);sin(theta)*sin(phi);cos(phi)];
e3 = -[cos(theta)*cos(phi);sin(theta)*cos(phi);-sin(phi)];
e2 = [-sin(theta);cos(theta);0];

b1 = p.b1;
b2 = p.b2;
b3 = p.b3;
% 
% a1 = lambda_1 + lambda_2 + lambda_3;
% a2 = -(lambda_1 * lambda_2 + lambda_2*lambda_3 + lambda_1*lambda_3);
% a3 = lambda_1 * lambda_2 * lambda_3;
% 
% b2 = 2.5; b3 =1;
B = [b1; b2; b3];
% e1'*B
T = [e1,e2,e3];
switch flag
    case 'Case1'
    lambda_1 = p.lambda_1;
   
    lambda_2 = p.lambda_2;
 
    lambda_3 = p.lambda_3;
    
    A = [lambda_1, 0 , 0; 0 , lambda_2, 0; 0 , 0 , lambda_3];

    case 'Case2'
    alpha = 0.5*(p.lambda_1 + p.lambda_2);
   
    beta = -0.5*(p.lambda_1 - p.lambda_2)*1j;
 
    lambda_3 = p.lambda_3;
    
    A = [alpha, beta, 0; -beta , alpha, 0; 0 , 0 , lambda_3];

end
A = T\A*T;

% C = e1'*T;
C = [1 0 0];
B = T\B;
R = eye(length(C)) - B*C*A;

% Jacobian of the sliding vector field
A_s = (eye(length(C))-(B*C*A)/(C*A*B))*A;
