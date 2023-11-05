function [A,B,C,phi2,theta2] = Jordan2OB(p)

% In OB form: three eigenvalues;
% B = [b1 b2 b3]
% C = [cos(theta)*sin(phi) sin(theta)*sin(phi) cos(phi)]
% are given with the C'T = 0

lambda_1 = p.lambda_1;

lambda_2 = p.lambda_2;

lambda_3 = p.lambda_3;

b1 = p.b1;

b2 =p.b2;

b3 = p.b3;

theta = p.theta;

phi = p.phi;



A = [lambda_1 0 0;0 lambda_2 0; 0 0 lambda_3];

B = [b1;b2;b3];
C = [cos(theta)*sin(phi) sin(theta)*sin(phi) cos(phi)];

% define the transformation explicitly with eigenspace


T = [1/(lambda_2*lambda_3),1/(lambda_1*lambda_3),1/(lambda_1*lambda_2);...
    -(lambda_2+lambda_3)/(lambda_2 *lambda_3),-(lambda_1+lambda_3)/(lambda_1*lambda_3),-(lambda_1+lambda_2)/(lambda_1* lambda_2);...
    1,1,1];

A = T*A/T;

B= T*B;

C = C/T;

phi2 = acos(C(3)/norm(C));
theta2 = atan(C(2)/C(1));
B = B*norm(C);
C = C/norm(C);






