% co-dimension set: {lambda_1, lambda_2, lambda_3, b2, b3, varphi, theta}

% 
lambda_1 = p.lambda_1;
lambda_2 = p.lambda_2;
lambda_3 = p.lambda_3;
b2       = p.b2;
b3       = p.b3;
varphi = p.varphi;
theta = p.theta;

%
a11 = -(cos(phi) - 1)*(cos(phi) + 1)*(lambda_1 - lambda_2)*cos(theta)^2 + (-lambda_2 + lambda_3)*cos(phi)^2 + lambda_2;
a12 = ((lambda_1 - lambda_2)*cos(theta)^2 + lambda_2 - lambda_3)*sin(phi)*cos(phi);
a13 = -cos(theta)*sin(theta)*sin(phi)*(lambda_1 - lambda_2);
a21 = ((lambda_1 - lambda_2)*cos(theta)^2 + lambda_2 - lambda_3)*sin(phi)*cos(phi);
a22 = ((lambda_1 - lambda_2)*cos(theta)^2 + lambda_2 - lambda_3)*cos(phi)^2 + lambda_3;
a23 = -cos(theta)*sin(theta)*cos(phi)*(lambda_1 - lambda_2);
a31 = -cos(theta)*sin(theta)*sin(phi)*(lambda_1 - lambda_2);
a32 = -cos(theta)*sin(theta)*cos(phi)*(lambda_1 - lambda_2);
a33 = (-lambda_1 + lambda_2)*cos(theta)^2 + lambda_1;
%
A = [a11 a12 a13; a21 a22 a23; a31 a32 a33];
%
C = [1,0,0];
B = [0;b2;b3];

R = eye(size(A,1)) - B*C*A;

