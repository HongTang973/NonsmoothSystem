function [A,B,C,phi2,theta2] = OB2Jordan(p)

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

a1 = lambda_1 + lambda_2 + lambda_3;
a2 = -(lambda_1 * lambda_2 + lambda_2*lambda_3 + lambda_1*lambda_3);
a3 = lambda_1 * lambda_2 * lambda_3;

A = [a1 1 0;a2 0 1; a3 0 0];

B = [b1;b2;b3];
C = [cos(theta)*sin(phi) sin(theta)*sin(phi) cos(phi)];

% define the transformation explicitly with eigenspace


T = [1/(lambda_2*lambda_3),1/(lambda_1*lambda_3),1/(lambda_1*lambda_2);...
    -(lambda_2+lambda_3)/(lambda_2 *lambda_3),-(lambda_1+lambda_3)/(lambda_1*lambda_3),-(lambda_1+lambda_2)/(lambda_1* lambda_2);...
    1,1,1];
%% do the transform 
[V,D]=eig(A);
eig_list = diag(D);
check = abs(eig_list.'-eig_list')>0;
n1  =0.5*length(check(check>0))
n2  =length(check(check==0))

% 
Tr = eye(n2);
P = [1 1;1i -1i];
for i=1:n1
    Tr = blkdiag(P,Tr);
    
end

% check 
Q=T*inv(Tr);
% inv(Q)*A*Q

% T = T*Tr;

A = real(Q\A*Q);

B= Q\B;

C = C*Q;

phi2 = real(acos(C(3)/norm(C)));
theta2 = real(atan(C(2)/C(1)));

e1 = [cos(theta2)*sin(phi2);sin(theta2)*sin(phi2);cos(phi2)];

B = B*norm(C);
C = C/norm(C);

if C*e1 <0
    theta2 = theta2 +pi;
end




