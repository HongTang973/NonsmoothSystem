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
lambda = -0.1;

A_1=[alpha beta 0; -beta alpha 0;0 0 lambda];
% A_1=[-3.2 -1 0; 25.61 0 -1;-75.03 0 0];

diag(eig(A_1))
[V,D]=eig(A_1)
theta = 0;
P = [cos(theta),0,-sin(theta);0,1,0;sin(theta),0,cos(theta)];
A_1 = inv(P)*A_1*P;
[V1,D1]=eig(A_1)
%
state0 = [-1 0 0]';
% 
f = -A_1*state0;
x0=0;
C=1;
y0=(C-lambda)/(beta*cos(theta));
z0=-1000;%-cot(theta)+(-exp(-pi*lambda/beta) + 1)/(2*sin(theta)*cos(theta));
InitCond=[x0 y0 z0]';
% InitCond=[0.142885021489076,0]';
% define the discontinuity manifold (plane)
% H(x) = 0
 C=[1 0 0];
% C*x=0;

% define the reset map
% when the event happens, impose the reset map on the plane (including 
% stretch and rotate)

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
T=20;
tspan=[0 T];
fs=100;
[tout,yout,teout,yeout,yeout0,ieout,amplitude]=PWSC_int(A_1,InitCond,state0,tspan,fs);
Period=teout(2:end)-teout(1:end-1);
% analytic results' solution
t=0:0.01:pi/beta;
x0
y0
z0
X=exp(lambda*t).*cos(beta*t) + exp(lambda*t).*sin(beta*t)*y0 - 1;
Y=-exp(lambda*t).*sin(beta*t) + exp(lambda*t).*cos(beta*t)*y0;
Z=exp(lambda*t)*z0;
figure
plot(t,X)
figure
plot(t,Y)
figure
plot(t,Z)

figure
plot(tout,yout(:,1))
if ~isempty(teout)
hold on
plot(teout,yeout(:,1),'b.')
end
figure
plot(tout,yout(:,2))
if ~isempty(teout)
hold on
plot(teout,yeout(:,2),'b.')
end
figure
plot(tout,yout(:,3))
if ~isempty(teout)
hold on
plot(teout,yeout(:,3),'b.')
end
% plot the orbit
S_T =50;
figure
plot3(yout(end-S_T*fs:end,1),yout(end-S_T*fs:end,2),yout(end-S_T*fs:end,3))
xlabel('x')
ylabel('y')
zlabel('z')
grid on
%
%
figure
plot3(yout(:,1),yout(:,2),yout(:,3))
%
figure
plot(yout(:,1),yout(:,2),'r-')
hold on
plot(yeout(:,1),yeout(:,2),'b.')



% script function
function [tout,yout,teout,yeout,yeout0,ieout,amplitude]=PWSC_int(A_1,InitCond,state0,tspan,fs)
%*************************#####************************* #
% checked and modified by Peter Tang /11th,3,2021       *#
%#************************#####************************* #
% gap:  backlash amount
% define integration timestep and timespan
tstart = tspan(1);
tfinal = tspan(2);
dt_scale = 1000/fs;

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
phi_1 = @(r_,z_) (1+4)*r_;
phi_2 = @(r_,z_) 0;
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
        if  r_( (ye+P_l3*r_(ye,a12,a13,f_(1))), a12, a13,f_(1)) < r_(ye,a12,a13,f_(1))
            % in the right direction to flip
        else % reverse the reset direction
            P_l3 = -P_l3;
        end
    end
    
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
    yeout = [yeout; ye];
    ieout = [ieout; ie];
    
   
%     A_1*ye'
    % Set the new initial conditions for next integration
    y0     = R(y(nt,:)',a12,a13,f_(1));
    tstart = t(nt)+1/fs/10;
    %
    
  

end

amplitude=1;

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
value = y(1);           %  detect switching point
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