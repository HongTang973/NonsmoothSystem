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
%
state0 = [-1 0 0]';
%
f = -A_1*state0;
x0=0;
C=2;
%  y0=(C-lambda)/(beta*cos(theta));
% z0=-(1-exp(-pi*lambda/beta) - 2*cos(theta)^2)/(2*sin(theta)*cos(theta))-1000;
% y0=2;
z0=10;
y0 = ((-alpha + lambda)*cos(theta)^2 + z0*(alpha - lambda)*sin(theta)*cos(theta) + C - lambda)/(beta*cos(theta));

InitCond=[x0 y0 z0]';
[1 0 0]*A_1*(InitCond-state0)
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

% %
% W=[0   0    0;
%    0   0    0;
%    0   0    0];

%
% B=W(:,1);
% A_s=(eye(3)-(B*C*A)/(C*A*B))*A;
% R = eye(3)+W;
% eig(A_s)
%
T=1000;
tspan=[0 T];
fs=100;
[tout,yout,teout,yeout,yeout0,ieout,amplitude]=PWSC_int(A_1,InitCond,state0,tspan,fs);
Period=beta*(teout(2:end)-teout(1:end-1))/pi;
teout = teout*0.2/pi;
tout = tout*0.2/pi;
% analytic results' solution
t=0:0.01:5*pi/beta;
% x0
% y0
% z0
X=cos(theta)^2*exp(lambda*t).*cos(beta*t) + sin(theta)^2*exp(lambda*t) + cos(theta)*exp(lambda*t).*sin(beta*t)*y0 + (cos(theta)*exp(lambda*t).*cos(beta*t)*sin(theta) - sin(theta).*exp(lambda*t)*cos(theta))*z0 - 1;
Y=-cos(theta)*exp(lambda*t).*sin(beta*t) + exp(lambda*t).*cos(beta*t)*y0 - sin(theta)*exp(lambda*t).*sin(beta*t)*z0;
Z=cos(theta)*exp(lambda*t).*cos(beta*t)*sin(theta) - sin(theta)*exp(lambda*t)*cos(theta) + sin(theta)*exp(lambda*t).*sin(beta*t)*y0 + (sin(theta)^2*exp(lambda*t).*cos(beta*t) + cos(theta)^2*exp(lambda*t))*z0;
% % X=0.2499999998*exp(-0.1*t) + 0.7500000004*exp(-0.1*t).*cos(0.2*t) + 0.8660254040*exp(-0.1*t).*sin(0.2*t)*y0 + (0.4330127018*exp(-0.1*t).*cos(0.2*t) - 0.4330127019*exp(-0.1*t))*z0 - 1;
% % Y=-0.8660254040*exp(-0.1*t).*sin(0.2*t) + (-4.*10^(-10)*exp(-0.1*t) + 1.000000000*exp(-0.1*t).*cos(0.2*t))*y0 - 0.5000000000*exp(-0.1*t).*sin(0.2*t)*z0;
% % Z=0.4330127020*exp(-0.1*t).*cos(0.2*t) - 0.4330127019*exp(-0.1*t) + 0.5000000000*exp(-0.1*t).*sin(0.2*t)*y0 + (0.7500000000*exp(-0.1*t) + 0.2500000000*exp(-0.1*t).*cos(0.2*t))*z0;
% figure
% plot(t,X)
% figure
% plot(t,Y)
% figure
% plot(t,Z)

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
% plot the orbit
S_T = 50;
figure(4)
plot3(yout(end-S_T*fs:end,1),yout(end-S_T*fs:end,2),yout(end-S_T*fs:end,3))
xlabel('x')
ylabel('y')
zlabel('z')
grid on
%
figure(5)
plot3(yout(:,1),yout(:,2),yout(:,3))
hold on
plot3(yeout(:,1),yeout(:,2),yeout(:,3),'r*')
a12=A_1(1,2);a13=A_1(1,3);
f_ = -A_1*state0;f1=f_(1);
V_l3 = -[0;a13;-a12]/sqrt(a12^2+a13^2);
P_l3 = [0;a12; a13]/sqrt(a12^2+a13^2);
r_ = @(y,a12,a13,f1) abs(a12*y(2)+a13*y(3)+f1)/sqrt(a12^2+a13^2);
tips = [max(yout(:,3))*V_l3';min(yout(:,3))*V_l3'];
plot3(tips(:,1),tips(:,2),tips(:,3),'r-')
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
 P=inv(V')*inv(V);
Check = A_1*P+P*A_1';
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


% Plot the contour line of the Lyapunov function
P_b = P(2:end,2:end);
V_L = 2;

% Lyp = @(x,y,z) P_b(1,1)*y.^2 + P_b(2,2)*z.^2 + 2*P_b(1,2)*y.^2 .*z.^2-V_L;
[x,y,z]=meshgrid(-10:0.2:10);
val=P_b(1,1)*y.^2 + P_b(2,2)*z.^2 + 2*P_b(1,2)*y.^2 .*z.^2;
xslice = [0];   
yslice = [];
zslice = [];
figure(4)
hold on
contourslice(x,y,z,val,xslice,yslice,zslice)
view(3)
grid on
% figure(4)
% x=0.001*x;
% isosurface(x,y,z,val,V_L)
% hold on
% figure(4)
% isosurface(x,y,z,val,0.5)
% figure(4)
% isosurface(x,y,z,val,1)
% fimplicit(@(y,z) P_b(1,1)*y.^2 + P_b(2,2)*z.^2 + 2*P_b(1,2)*y.^2 .*z.^2-V_L,[-10 10 -10 10])
% ezimplot3(Lyp,[-10,10]);
% Lyp(0,1,2)
% figure
% plot(tout,Deriv)
% figure
% plot(tout,K_Energy)
% figure
% plot(tout,P_Energy)
% figure
% plot(tout,K_Energy+P_Energy)
% script function
function [tout,yout,teout,yeout,yeout0,ieout,amplitude]=PWSC_int(A_1,InitCond,state0,tspan,fs)
%*************************#####************************* #
% checked and modified by Peter Tang /11th,3,2021       *#
%#************************#####************************* #
% gap:  backlash amount
% define integration timestep and timespan
tstart = tspan(1);
tfinal = tspan(2);
% dt_scale = 1000/fs;
amplitude=[];
% define initial condition to start
y0 = InitCond;

%define options
refine = 1;
options = odeset('Events',@(t,y) DetectingZero_events(t,y),'RelTol',1e-12,'AbsTol',1e-12,...
    'Refine',refine);

% start calculation
tout = tstart;
yout = y0.';
teout = [];
yeout = [];
yeout0= [];
ieout = [];
f_ = -A_1*state0;


% define the reset map
a12=A_1(1,2);a13=A_1(1,3);
V_l3 = [0;a13;-a12]/sqrt(a12^2+a13^2);
P_l3 = [0;a12; a13]/sqrt(a12^2+a13^2);
r_ = @(y,a12,a13,f1) abs(a12*y(2)+a13*y(3)+f1)/sqrt(a12^2+a13^2);
z_ = @(y,V_l3) V_l3'*y;
%
phi_1 = @(r_,z_) (1+0.8)*r_;
phi_2 = @(r_,z_) 8*z_;
%
R = @(y,a12,a13,f1) y + phi_1(r_(y,a12,a13,f1),z_(y,V_l3))*P_l3 + phi_2(r_(y,a12,a13,f1),z_(y,V_l3))*V_l3;

while tout(end)<tfinal
    % define out put time serias
    timespan=[tstart:1/fs:tfinal];
    %     timespan=[tstart tfinal];
    
    % exit case 1
    if length(timespan)<2
        % last time step
        break;
    end
    % change the stiffness matrix with different region
    A_matrix=A_1;
    %         offset=0;
    % Solve until the first terminal event.
    [t,y,te,ye,ie] = ode45(@(t,y) f(t,y,A_matrix,state0),...
        timespan,y0,options);
    
    % correct the right direction
    if ~isempty(ye)&& a12*ye(2)+a13*ye(3)+f_(1)<0 % Incoming point
        try
            if  r_( (ye+P_l3*r_(ye,a12,a13,f_(1))), a12, a13,f_(1)) < r_(ye,a12,a13,f_(1))
                % in the right direction to flip
            else % reverse the reset direction
                P_l3 = -P_l3;
            end
        catch ME
            figure
            plot(tout,yout(:,1))
            figure
            plot(tout,yout(:,2))
            figure
            plot(tout,yout(:,3))
            keyboard
        end
    end
    r1 = r_(y0,a12,a13,f_(1));
    z1 = z_(y0,V_l3);
    % Accumulate output.  This could be passed out as output arguments.
    nt    = length(t);
    
    % exit case II
    %     if max(abs(yout(:,1)))>10 % divergence or failure
    %         % last time step
    %         break;
    %     end
    
    % output history
    tout  = [tout; t(2:nt)];
    yout  = [yout; y(2:nt,:)];
    % output events
    teout = [teout; te];
    yeout0 = [yeout0; y0'];
    yeout = [yeout; ye];
    ieout = [ieout; ie];
    
    if ~isempty(ye)
        r2 = r_(ye,a12,a13,f_(1));
         z2 = z_(ye',V_l3);
         try
        amplitude=[amplitude,[r2/r1;z2/z1]];
         catch ME
             keyboard
         end
    end
    
    %     A_1*ye'
    % Set the new initial conditions for next integration
    A_1*(ye'-state0)
    y0     = R(y(nt,:)',a12,a13,f_(1));
    y0(1)  = 0;
    tstart = t(nt)+1/fs/10;
    %
    
    
    
end



%------------------------------------------------------------------------%
% function during region 1
    function dydt = f(t,y,A_matrix,state0)
        
        dydt = A_matrix*(y-state0);
        
    end

%------------------------------------------------------------------------%
% function reset map


%------------------------------------------------------------------------%
    function [value,isterminal,direction] = DetectingZero_events(t,y)
        % Locate the time when flap hit through boundary in both directions
        % and stop integration when isterminal ==1
        % Events: switch boundary; zero heave velocity; zero pitch velocity; zero
        % flap velocity
        value = (y(1));           %  detect switching point
        isterminal = 1;                           %  stop the integration
        direction  = -1;                           %  mark with event function
        %  adding direction
    end
%
end

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