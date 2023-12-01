% to construct the 3d hybrid system
% the reset map will expand the Lyapunov function level set
% Case I: equilibrium on the boundary
%        the eigenvector v3 is on the surface H(x) = 0
% generally we can always rotate the eigenvector of the real one to be in
% the plane XoZ; specifically, in this case relative angle theta to the Z axis.
clc
clear
close all
% define the linear part
alpha = -0.1;
beta  = 0.2;
lambda = -0.5;

A_1=[alpha beta 0; -beta alpha 0;0 0 lambda];

diag(eig(A_1))
[V,D]=eig(A_1)
theta = pi/6;
P = [cos(theta),0,-sin(theta);0,1,0;sin(theta),0,cos(theta)];
A_1 = inv(P)*A_1*P;
[V1,D1]=eig(A_1)
%  transfrom into the real coordinate 
RV = [];
RV(:,1) = V1(:,1);
RV(:,2) = 0.5*(V1(:,2)+V1(:,3));
RV(:,3) = -0.5i*(V1(:,2)-V1(:,3))
    
%
state0 = [-1 0 0]';
%
f = -A_1*state0;
x0=0;
C=0.1;
%  y0=(C-lambda)/(beta*cos(theta));
% z0=-(1-exp(-pi*lambda/beta) - 2*cos(theta)^2)/(2*sin(theta)*cos(theta))-1000;
% y0=2;
z0=2;
y0 = ((-alpha + lambda)*cos(theta)^2 + z0*(alpha - lambda)*sin(theta)*cos(theta) + C - lambda)/(beta*cos(theta));
%  y0 =1.22;z0=0;
% InitCond=[x0 y0 z0]';
% 
InitCond = [0;15.5437868298994;11.3200980979380];
% InitCond=[0;1.652370333976446;0.130193205604953 ];
vx=[1 0 0]*A_1*(InitCond-state0)
% Boundary_equilirium = -inv(A_1(2:3,2:3))*f(2:3)
% Special case: initial condition on the eigenvector v3
% InitCond=state0+5*V1(:,3);
% InitCond=[0,0.5,2.5205]';
% InitCond=[0.142885021489076,0]';
% define the discontinuity manifold (plane)
% H(x) = 0
C=[1 0 0];
% C*x=0;
n_p = [1,0,0]';c_p=0;n_L = V1(:,3);origin_L=state0;
pos = get_intersection_point(n_p,c_p,n_L,origin_L);
% InitCond=pos+[0,0.01,];
% define the reset map
% when the event happens, impose the reset map on the plane (including stretch and rotate)
pr=0.8; %coef of restitution
pz = 4.6;
%
gamma = (alpha-lambda)/beta;
eta = lambda/beta;
psy =eta*(-gamma*(1 + pr)*cos(theta)^2 + pz*(eta*gamma + gamma^2 - 1)*sin(theta) ...
    + (1 + pr)*(eta + 2*gamma))/((gamma*(1 + pr)*(eta*gamma + gamma^2 + 1)*cos(theta)^2 ...
    - (gamma^2 + 1)*(-sin(theta)*pz + (1 + pr)*(gamma + eta)))*cos(theta))

psz=((-gamma*(1 + pr)*(eta*gamma + gamma^2 + 1)*sin(theta) + ...
    pz*(gamma^2 + 1))*cos(theta)^2 - eta*((1 + pr)*(eta*gamma + gamma^2 - 1)*sin(theta) ...
    - pz*(eta + 2*gamma)))/((gamma*(1 + pr)*(eta*gamma + gamma^2 + 1)*cos(theta)^2 ...
    - (gamma^2 + 1)*(-sin(theta)*pz + (1 + pr)*(gamma + eta)))*cos(theta))

% %
(gamma*(1 + pr)*(eta*gamma + gamma^2 + 1)*sin(theta)^2 + (-gamma^2*pz - pz)*sin(theta) + eta*(1 + pr))*beta/((gamma^2*sin(theta)^2 + 1)*(1 + pr))
(gamma*(1 + pr)*(eta*gamma + gamma^2 + 1)*cos(theta)^2 - (-sin(theta)*pz + (1 + pr)*(gamma + eta))*(gamma^2 + 1))*cos(theta)/((gamma*(1 + pr)*(eta*gamma + gamma^2 + 1)*sin(theta) - gamma^2*pz - pz)*cos(theta)^2 + eta*((1 + pr)*(eta*gamma + gamma^2 - 1)*sin(theta) - pz*(eta + 2*gamma)))
% W=[0   0    0;
%    0   0    0;
%    0   0    0];
pz0=((1 + pr)*((gamma + eta)*(gamma^2*sin(theta)^2 + 1) - gamma*cos(theta)^2))/((gamma^2 + 1)*sin(theta))
%
% B=W(:,1);
% A_s=(eye(3)-(B*C*A)/(C*A*B))*A;
% R = eye(3)+W;
% eig(A_s)
%
T=800;
tspan=[0 T];
fs=100;
mu =1;
[tout,yout,teout,yeout,yeout0,ieout,amplitude]=PWSC_int(A_1,InitCond,state0,tspan,fs,mu,0);
Period=beta*(teout(2:end)-teout(1:end-1))/pi;
figure(100)
plot(tout,yout(:,1))
hold on
figure(99)
plot(yout(end-50*fs:end,1),yout(end-50*fs:end,3))
% use the analytical method 
ana_yout=[];
for tt=1:length(tout)
    DS=diag(exp(diag(D1)*tout(tt)));
    y_temp = V1*DS*inv(V1)*(InitCond-state0)+state0;
    ana_yout=[ana_yout;y_temp'];
end
figure(100)
plot(tout,ana_yout(:,1))
T = Period(end)*pi/beta;
% construct the Jacobi matrix of the composed map
a11=A_1(1,1);a12=A_1(1,2);a13=A_1(1,3);
f_ = -A_1*state0;f1=f_(1);
N = sqrt(a12^2+a13^2);
V_l3 = [0;a13;-a12]/N;
P_l3 = [0;a12; a13]/N;
W =  -(1+pr)*P_l3*[a11 a12 a13]/N+pz*V_l3*[a11 a12 a13]/N;
 
T_2=0;
T_1=T-T_2;
f_ls= @(t,y) A_1*(y-state0);
DS=diag(exp(diag(D1)*T_1));
DS0=diag(exp(diag(D1)*T_2));
S=diag(DS);
S0=diag(DS0);
KK0=V1*DS0*inv(V1);
KK1=V1*DS*inv(V1);
x_00 = yeout0(end,:)';
x_0  = yeout(end,:)';
% x_00 = [0 15.5398,11.3170]';
% x_0 = [0 -3.2528,-0.5722]';
x_00(1)=0;
x_0(1)=0;
C =[1 0 0];
P = eye(3) + W;
x_00 = P*(x_0-state0)+state0
correction = ((f_ls(0,x_00)-P*f_ls(0,x_0))*C)/(C*f_ls(0,x_0));
PP=P+1*correction;
J=PP*KK1;
[JV,JD]=eig(J)
JV(:,1)'*(x_00-state0)/norm(JV(:,1))/norm(x_00-state0)
JV(:,2)'*(x_00-state0)/norm(JV(:,2))/norm(x_00-state0)
teout = teout*0.2/pi;
tout = tout*0.2/pi;
% analytic results' solution
t=0:0.01:5*pi/beta;
% x0
% y0
% z0
% X=cos(theta)^2*exp(lambda*t).*cos(beta*t) + sin(theta)^2*exp(lambda*t) + cos(theta)*exp(lambda*t).*sin(beta*t)*y0 + (cos(theta)*exp(lambda*t).*cos(beta*t)*sin(theta) - sin(theta).*exp(lambda*t)*cos(theta))*z0 - 1;
% Y=-cos(theta)*exp(lambda*t).*sin(beta*t) + exp(lambda*t).*cos(beta*t)*y0 - sin(theta)*exp(lambda*t).*sin(beta*t)*z0;
% Z=cos(theta)*exp(lambda*t).*cos(beta*t)*sin(theta) - sin(theta)*exp(lambda*t)*cos(theta) + sin(theta)*exp(lambda*t).*sin(beta*t)*y0 + (sin(theta)^2*exp(lambda*t).*cos(beta*t) + cos(theta)^2*exp(lambda*t))*z0;
% % X=0.2499999998*exp(-0.1*t) + 0.7500000004*exp(-0.1*t).*cos(0.2*t) + 0.8660254040*exp(-0.1*t).*sin(0.2*t)*y0 + (0.4330127018*exp(-0.1*t).*cos(0.2*t) - 0.4330127019*exp(-0.1*t))*z0 - 1;
% % Y=-0.8660254040*exp(-0.1*t).*sin(0.2*t) + (-4.*10^(-10)*exp(-0.1*t) + 1.000000000*exp(-0.1*t).*cos(0.2*t))*y0 - 0.5000000000*exp(-0.1*t).*sin(0.2*t)*z0;
% % Z=0.4330127020*exp(-0.1*t).*cos(0.2*t) - 0.4330127019*exp(-0.1*t) + 0.5000000000*exp(-0.1*t).*sin(0.2*t)*y0 + (0.7500000000*exp(-0.1*t) + 0.2500000000*exp(-0.1*t).*cos(0.2*t))*z0;
% figure
% plot(t,X)
% figure
% plot(t,Y)
% figure
% plot(t,Z)
S_T = 499;
figure(1)
plot(tout,yout(:,1))
hold on
% plot(t,X,'k--','linewidth',1.2)
if ~isempty(yeout)
    hold on
    plot(teout,yeout(:,1),'b.')
    
end
figure(2)
plot(tout,yout(:,2))
hold on
% plot(t,Y,'k--','linewidth',1.2)
if ~isempty(yeout)
    hold on
    plot(teout,yeout(:,2),'b.')
end
figure(3)
plot(tout,yout(:,3))
hold on
% plot(t,Z,'k--','linewidth',1.2)
if ~isempty(yeout)
    hold on
    plot(teout,yeout(:,3),'b.')
end
%
% Chose a positive definite matrix Q to find corresponding P matrix to
% construct the Lyapunov function
[V,D]=eig(A_1);
Q=diag(abs(real(diag(D))));
Q=diag([1 2 4 ]);
P=lyap(A_1',Q);
% P=inv(V')*inv(V);
% Check = A_1*P+P*A_1';
%
% Plot the contour line of the Lyapunov function
P_b = real(P);
V_L = 2;

% Lyp = @(x,y,z) P_b(1,1)*y.^2 + P_b(2,2)*z.^2 + 2*P_b(1,2)*y.^2 .*z.^2-V_L;
x_0 =1;
[x,y,z]=meshgrid(-10:0.2:10);
val=P_b(1,1)*(x+x_0).^2 + P_b(2,2)*y.^2 + P_b(3,3)*z.^2 + 2*P_b(1,2)*y.*(x+x_0)+ 2*P_b(1,3)*z.*(x+x_0)+ 2*P_b(2,3)*y.*z;
xslice = [0];   
yslice = [];
zslice = [];
lvls = 0:10:100;

% Implicitplot

Imf = @(x,y,z) P_b(1,1)*(x+x_0).^2 + P_b(2,2)*y.^2 + P_b(3,3)*z.^2 + 2*P_b(1,2)*y.*(x+x_0)+ 2*P_b(1,3)*z.*(x+x_0)+ 2*P_b(2,3)*y.*z -1;
interval = [-10 10 -10 10 -10 10];
figure
fimplicit3(Imf,interval)
xlabel('x')
ylabel('y')
zlabel('z')


% plot the orbit

figure(4)
plot3(yout(end-S_T*fs:end,1),yout(end-S_T*fs:end,2),yout(end-S_T*fs:end,3))
xlabel('x')
ylabel('y')
zlabel('z')
grid on
hold on
contourslice(x,y,z,val,xslice,yslice,zslice,[5,lvls])
view(3)
grid on
%
figure(5)
plot3(yout(:,1),yout(:,2),yout(:,3))
hold on
plot3(yeout(:,1),yeout(:,2),yeout(:,3),'r*')

tips = [max(yout(:,3))*V_l3';min(yout(:,3))*V_l3'];
plot3(tips(:,1),tips(:,2),tips(:,3),'r-')
grid on
contourslice(x,y,z,val,xslice,yslice,zslice,[5,lvls])
view(3)
% for j=1:length(yeout(:,1))
%     figure(5)
%     point = yeout(j,:)+r_(yeout(j,:),a12,a13,f(1))*P_l3';
%     tips = [yeout(j,:);point];
%     plot3(tips(:,1),tips(:,2),tips(:,3),'r-')
% end
%
figure(6)
plot(yout(:,1),yout(:,2),'r-')
if ~isempty(yeout)
    hold on
    plot(yeout(:,1),yeout(:,2),'b.')
end


%
figure(7)
plot(yout(end-S_T*fs:end,1),yout(end-S_T*fs:end,3),'r-')
if ~isempty(yeout)
    hold on
   % plot(yeout(:,1),yeout(:,2),'b.')
end
% Chose a positive definite matrix Q to find corresponding P matrix to
% construct the Lyapunov function
[V,D]=eig(A_1);
Q=diag(abs(real(diag(D))));
Q=diag([1 2 3 ]);
% P=lyap(A_1',Q);
%  P=inv(V')*inv(V);
% Check = A_1*P+P*A_1';
%

%% inspect the constructed Lyapunov function along the trajectory
focus=state0;
Yout=yout(end-S_T*fs:end,:);
Lyp=zeros(size(Yout,1),1);
%Deriv=zeros(size(yout,1),1);
%K_Energy=zeros(size(yout,1),1);
%P_Energy=zeros(size(yout,1),1);
for kk= 1:size(Yout,1)
    Yout(kk,:)=Yout(kk,:)-focus';
    Lyp(kk)=Yout(kk,:)*P*Yout(kk,:)';
    %Deriv(kk)=-Yout(kk,:)*Q*Yout(kk,:)';
    %K_Energy(kk)= 0.5*yout(kk,4:6)*M_new*yout(kk,4:6)'+0.*yout(kk,1:3)*K_new*yout(kk,1:3)';
    %P_Energy(kk)= 0.5*(yout(kk,1:3)-focus(1:3)')*K_new*(yout(kk,1:3)-focus(1:3)')'-0*0.5*focus(1:3)'*K_new*focus(1:3);
end
figure(8)
plot(tout(end-S_T*fs:end,1),Lyp)
xlabel('t')
title('Lyapunov function')

% plot the time series in the eigenvet coordinate

cyout = P*yout';
cyeout = P*yeout';

figure
plot(tout,cyout(1,:))
hold on
plot(teout,cyeout(1,:),'b.')
figure
plot(tout,cyout(2,:))
hold on
plot(teout,cyeout(2,:),'b.')
figure
plot(tout,cyout(3,:))
hold on
plot(teout,cyeout(3,:),'b.')

% get the intersection point between a line and a plane(3d)
function pos = get_intersection_point(n_p,c_p,n_L,origin_L)
% the parameter equation of the line
% pos = origin_L + n_L*t;
% the equation of the plane:
% n_p' *  coord = c_p
% so
t = (c_p -n_p'*origin_L)/(n_p'*n_L);
pos = origin_L + n_L*t;
end