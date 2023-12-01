% Numerically calculate the discrepancy of the Lyapunov function 
% Independent variable: the state variable, on the discontinuity surface
% of the system 
% PS: from the positive set
% LF_1: the decrease in the flow region
% LF_2: the elevaltion by the reset map(changing the pr and pz)
%*************** start parallel computation ***************%
% find the numbers of this 
Initialization_for_Par;

% define the system
clc
clear
close all
% define the linear part
alpha = -0.1;
beta  = 0.2;
lambda = -0.5;
gamma = (alpha-lambda)/beta;
A_1=[alpha beta 0; -beta alpha 0;0 0 lambda];
theta = pi/6;
P = [cos(theta),0,-sin(theta);0,1,0;sin(theta),0,cos(theta)];
A_1 = inv(P)*A_1*P;

% define initial condition: grid on the positive set
%
state0 = [-1 0 0]';
F_ls= @(t,y) A_1*(y-state0);
%
f = -A_1*state0;
x0=0;
% C=0.1:0.1:10;
% z0 = -20:0.1:20;
[C,Z]=meshgrid([[0.01:0.01:0.1],[0.1:0.1:10]],[-20:0.1:20]);
% on the parallel line with distance C
Y = ((-alpha + lambda)*cos(theta)^2 + Z*(alpha - lambda)*sin(theta)*cos(theta) + C - lambda)/(beta*cos(theta));
fs =100;
tspan = [0 100];
% define Lyapunov function
% Chose a positive definite matrix Q to find corresponding P matrix to
% construct the Lyapunov function
[V,D]=eig(A_1);
Q=diag([1 1 1 ]);
P_b=real(lyap(A_1',Q));
x_0 =1;
Lf = @(x,y,z) P_b(1,1)*(x+x_0).^2 + P_b(2,2)*y.^2 + P_b(3,3)*z.^2 + 2*P_b(1,2)*y.*(x+x_0)+ 2*P_b(1,3)*z.*(x+x_0)+ 2*P_b(2,3)*y.*z;
% define the reset map
a11=A_1(1,1);a12=A_1(1,2);a13=A_1(1,3);
N = sqrt(a12^2+a13^2);
V_l3 = [0;a13;-a12]/N;
P_l3 = [0;a12; a13]/N;
pr=0.8; %coef of restitution
pz = 8;
vx_= @(t,y)  [1 0 0]*A_1*(y-state0); 
Wx_= (-(1+pr)*[0 a12 a13]'+pz*[0 a13 -a12]')./(a12^2+a13^2);
MW =  -(1+pr)*P_l3*[a11 a12 a13]/N+pz*V_l3*[a11 a12 a13]/N;
%Pseudo equilibrium
%
gamma = (alpha-lambda)/beta;
eta = lambda/beta;
psy =eta*(-gamma*(1 + pr)*cos(theta)^2 + pz*(eta*gamma + gamma^2 - 1)*sin(theta) ...
    + (1 + pr)*(eta + 2*gamma))/((gamma*(1 + pr)*(eta*gamma + gamma^2 + 1)*cos(theta)^2 ...
    - (gamma^2 + 1)*(-sin(theta)*pz + (1 + pr)*(gamma + eta)))*cos(theta));

psz=((-gamma*(1 + pr)*(eta*gamma + gamma^2 + 1)*sin(theta) + ...
    pz*(gamma^2 + 1))*cos(theta)^2 - eta*((1 + pr)*(eta*gamma + gamma^2 - 1)*sin(theta) ...
    - pz*(eta + 2*gamma)))/((gamma*(1 + pr)*(eta*gamma + gamma^2 + 1)*cos(theta)^2 ...
    - (gamma^2 + 1)*(-sin(theta)*pz + (1 + pr)*(gamma + eta)))*cos(theta));

%
stop_flag =1;
Dis = [];
ET =[];
Norm_eig = [];
tic
parfor i=1:size(Y,1)*size(Z,2)
    row = floor(i/size(Y,2))+1;
    col = mod(i,size(Y,2));
    if col ==0
        col =size(Z,2);
        row = row-1;
    end
    %
    InitCond=[x0,Y(row,col),Z(row,col)]';
    [~,~,teout,yeout,yeout0,~,~]=PWSC_int(A_1,InitCond,state0,tspan,fs,1,stop_flag);
    LF_1 = Lf(x0,Y(row,col),Z(row,col)) - Lf(yeout(1),yeout(2),yeout(3));
    % reset map
    Back = yeout'+ Wx_*vx_(0,yeout');
    LF_2 = Lf(Back(1),Back(2),Back(3)) - Lf(yeout(1),yeout(2),yeout(3));
    Dis(i) = LF_1 -LF_2;
    ET(i) = teout;
    Norm_eig(i)=FloqueMultiplier(A_1,teout,InitCond,yeout',F_ls,MW);
end
toc
%
LD      =zeros(size(C,1),size(C,2));
EVT     =zeros(size(C,1),size(C,2));
MainEig =zeros(size(C,1),size(C,2));
for i=1:size(C,1)*size(C,2)
    row = floor(i/size(Y,2))+1;
    col = mod(i,size(Y,2));
    if col ==0
        col =size(Z,2);
        row = row-1;
    end
    LD(row,col) = Dis(i);
    EVT(row,col) = ET(i);
    MainEig(row,col) = Norm_eig(i);
end
Z0 =[-20 20];
Y0 = ((-alpha + lambda)*cos(theta)^2 + Z0*(alpha - lambda)*sin(theta)*cos(theta) - lambda)/(beta*cos(theta));
%
save(sprintf('20220210-111%3fAlpha%3fBeta%3fLambda%3fpr%3fpz.mat',-alpha,beta,-lambda,pr,pz))
%
% figure;plot(Y0,Z0,'b')
figure(1)
hold on
plot(15.5398,11.3170,'ro','displayname','LCO-x0')
plot(-3.2528,-0.5722,'k*','displayname','LCO-x1')
plot(Y0,Z0,'b','displayname','zero velocity line')
plot(psy,psz,'gs','displayname','Pseudo Equilibrium')
%
(15.5398+3.2528)/(11.3170+0.5722)
-(sin(theta)*pz*gamma + pr + 1)/((1 + pr)*sin(theta)*gamma - pz)

xlabel('y')
ylabel('z')
legend('location','best')
figure(1)
points1=contour(Y,Z,LD,[0 0]);
points2=contour(Y,Z,MainEig,[0.9 1 1.2 1.4 1.6 1.8 1.86],'showtext','on');
 points3=contour(Y,Z,EVT,[16 15 10 5],'showtext','on');
function t=FloqueMultiplier(A_1,T,x_00,x_0,F_ls,MW)
%%  construct the Jacobi matrix of the composed map
W = eye(size(A_1,1)) + MW;
[V1,D1]=eig(A_1);
T_2=0;
T_1=T-T_2;

DS=diag(exp(diag(D1)*T_1));
DS0=diag(exp(diag(D1)*T_2));
% S=diag(DS);
% S0=diag(DS0);
KK0=V1*DS0*inv(V1);
KK1=V1*DS*inv(V1);
% x_00 = yeout0(end,:)';
% x_0  = yeout(end,:)';
C =[1 0 0];
correction = ((F_ls(0,x_00)-W*F_ls(0,x_0))*C)/(C*F_ls(0,x_0));
PP=W+0*correction;
J=KK0*PP*KK1;
[~,JD]=eig(J);
t = max(abs(diag(JD)));
end