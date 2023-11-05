% co-dimension set: {alpha, lambda_2, lambda_3, b2, b3, varphi, theta}

% 
alpha= p.alpha;
beta = p.beta;
lambda_3 = p.lambda_3;
b2       = p.b2;
b3       = p.b3;
varphi = p.varphi;
theta = p.theta;

%
a11 = (-alpha + lambda)*cos(phi)^2 + alpha;
a12 = sin(phi)*cos(phi)*(alpha - lambda);
a13 = beta*sin(phi);
a21 = sin(phi)*cos(phi)*(alpha - lambda);
a22 = (alpha - lambda)*cos(phi)^2 + lambda;
a23 = beta*cos(phi);
a31 = -beta*sin(phi);
a32 = -beta*cos(phi);
a33 = alpha;
%
A = [a11 a12 a13; a21 a22 a23; a31 a32 a33];
%
C = [1,0,0];
B = [0;b2;b3];

R = eye(size(A,1)) - B*C*A;



