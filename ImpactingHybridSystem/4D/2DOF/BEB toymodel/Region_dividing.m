% define the linear part
alpha = -0.1;
beta  = 0.2;
lambda = -0.5;

A_1=[alpha beta 0; -beta alpha 0;0 0 lambda];
theta = pi/6;
P = [cos(theta),0,-sin(theta);0,1,0;sin(theta),0,cos(theta)];
A_1 = inv(P)*A_1*P;

% define initial condition: grid on the positive set
%
state0 = [-1 0 0]';

%
f = -A_1*state0;
% construct the Jacobi matrix of the composed map
a12=A_1(1,2);a13=A_1(1,3);
N = a12^2+a13^3;
W = eye(3) - (1+pr)*[0 0 0; 0 a12^2 a12*a13;0 a12*a13 a13^2]/N+pz*[0 0 0;0 a12*a13 a13^2;0 -a12^2 -a12*a13]/N;
% Chose a positive definite matrix Q to find corresponding P matrix to
% construct the Lyapunov function
[V,D]=eig(A_1);
Q=diag([1 1 1 ]);
Omega=real(lyap(A_1',Q));
P_b = -W'*(Omega'+Omega);

x_0 =1;
Direction = @(y,z) P_b(1,1)*(0+x_0).^2 + P_b(2,2)*y.^2 + P_b(3,3)*z.^2 + 2*P_b(1,2)*y.*(0+x_0)+ 2*P_b(1,3)*z.*(0+x_0)+ 2*P_b(2,3)*y.*z;
[C,Z]=meshgrid([[0.01:0.01:0.1],[-10:0.1:10]],[-20:0.2:20]);
% on the parallel line with distance C
Y = ((-alpha + lambda)*cos(theta)^2 + Z*(alpha - lambda)*sin(theta)*cos(theta) + C - lambda)/(beta*cos(theta));
Dire = Direction(Y,Z);
figure(1)
points=contour(Y,Z,Dire,[0 -500],'ShowText','on');
hold on