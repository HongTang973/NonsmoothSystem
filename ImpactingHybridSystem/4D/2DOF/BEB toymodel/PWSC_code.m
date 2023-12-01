clc
clear
close all
% to find the simplest BEB case in 3 or higher dimension system
% Example 5.8
% A_1=[-0.8 1 0;-0.57 0 1;-0.09 0 0];
% A_2=[-0.1 1 0;-0.2 0 1;-60 0 0];
% LCO example
A_1 =[-5 1 0;-9 0 1 ;-5 0 0];
% A_2 =A_1;
A_2 =[-5 1 0;12 0 1; -14 0 0];
% A_1 = A_2;
% manifold of 3-d system
% A_1=[-0.2 10 0;-10 -0.2 0; 0 0 -0.3];
% A_2=A_1;
% invariant cones
% A_1=[-1 -1 0;1.28 0 -1;-0.624 0 0];
% A_2=[-3.2 -1 0; 25.61 0 -1;-75.03 0 0];
[V1,ind1]=eig(A_1)
Q1=Jordan_T(A_1);
%  inv(Q1)*A_1*Q1
P1=inv(V1')*inv(V1)
%
[V2,ind2]=eig(A_2)
P2=inv(V2')*inv(V2)
%
Q2=Jordan_T(A_2);
inv(Q2)*A_2*Q2
vect_preload=[0 0 1]';
eq1=-A_1\vect_preload
eq2=-A_2\vect_preload
%
OB1=A_2'*P1+P1*A_2
eig(OB1)
%
OB2=A_1'*P2+P2*A_1
eig(OB2)
InitCond=0.0001*[1 1 1]';
%
% InitCond=eq2+0.01;
%
T=300;
%
tspan=[0 T];
%
fs=50;
[tout,yout,teout,yeout,yeout0,ieout,amplitude]=PWSC_int(A_1,A_2,InitCond,vect_preload,tspan,fs);
%
eig_yeout1=inv(Q1)*yeout';
eig_yeout2=inv(Q2)*yeout';
%
figure
plot(tout,yout(:,1))

figure
plot(tout,yout(:,2))

figure
plot(tout,yout(:,3))

figure
plot3(yout(:,1),yout(:,2),yout(:,3))
hold on
%
scale=10;
data_set = scale*[V1(:,3),-V1(:,3),V2(:,3),-V2(:,3)];
X=reshape(data_set(1,:),2,2);
Y=reshape(data_set(2,:),2,2);
Z=reshape(data_set(3,:),2,2);
plot3(X,Y,Z,'o-');
plot3(eq1(1),eq1(2),eq1(3),'s')
plot3(yeout(end,1),yeout(end,2),yeout(end,3),'*')
plot3(yeout(end-1,1),yeout(end-1,2),yeout(end-1,3),'*')
xlabel('x')
ylabel('y')
zlabel('z')
%%
function Q=Jordan_T(A)
[V,D]=eig(A);
eig_list = diag(D);
check = abs(eig_list.'-eig_list')>0;
n1  =0.5*length(check(check>0));
n2  =length(check(check==0));
% remember the location of the real eigenvalue
ind=(check==0);
%
% 
Tr = [];
P = [1 1;1i -1i];
i=1;
while i<=(2*n1+n2)
    if ind(i)
    Tr = blkdiag(Tr,eye(1));
    i=i+1;
    else
    Tr = blkdiag(Tr,P);
    i=i+2;    
    end
end
% check 
Q=V*inv(Tr);
end
function [tout,yout,teout,yeout,yeout0,ieout,amplitude]=PWSC_int(A_1,A_2,InitCond,vect_preload,tspan,fs)
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
    if (y0(1))<0      %case 1
        A_matrix=A_1;
%         offset=0;
        options = odeset(options,'MaxStep',1e-3);
    elseif (y0(1))==0 % case 2: Exactly hit the boundary
        %
        if y0(2)   %Get into region2 in the next time step
            A_matrix=A_2;
%             offset=sign(y0(1));
            %         options = odeset(options,'InitialStep',t(end)-t(end-1));
            options = odeset(options,'InitialStep',dt_scale*1e-5,'MaxStep',dt_scale*1e-4);

        else               %Get into region1 in the next time step
            A_matrix=A_1;
%             offset=0;
            options = odeset(options,'MaxStep',1e-3);
        end
        %
    else %case3: Located at region 2
        A_matrix=A_2;
%         offset=sign(y0(1));
        options = odeset(options,'InitialStep',dt_scale*1e-5,'MaxStep',dt_scale*1e-4);
    end
    
    % Solve until the first terminal event.
    [t,y,te,ye,ie] = ode45(@(t,y) f(t,y,A_matrix,vect_preload),...
        timespan,y0,options);
    
    
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
    
   
   
    % Set the new initial conditions for next integration
    y0     = y(nt,:);
    yeout0=[yeout0;y0];
    tstart = t(nt);
    % A good guess of a valid first timestep is the length of the last valid
    % timestep, so use it for faster computation.  'refine' is 4 by default.
    %options = odeset(options,'InitialStep',t(nt)-t(nt-refine),...
    %  'MaxStep',t(nt)-t(1));

end

amplitude=1;

%------------------------------------------------------------------------%
% function during region 1
function dydt = f(t,y,A_matrix,vect_preload)

dydt = A_matrix*y+vect_preload;

end
%------------------------------------------------------------------------%
function [value,isterminal,direction] = DetectingZero_events(t,y)
% Locate the time when flap hit through boundary in both directions
% and stop integration when isterminal ==1
% Events: switch boundary; zero heave velocity; zero pitch velocity; zero
% flap velocity
value = [(y(1));abs(y(2))];           %  detect switching point
isterminal = [1;0];                           %  stop the integration
direction  = [0;0];                           %  mark with event function
                                              %  adding direction
end
%

end