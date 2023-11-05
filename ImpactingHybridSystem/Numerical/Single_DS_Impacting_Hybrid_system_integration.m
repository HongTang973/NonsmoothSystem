%  Numerical integration of the Impacting Hybrid system with:
% the A matrix: F_ls = A*y  defing the flow of the linear differential
% equation 
% R : the reset map matrix, this means the flow is continuous 
% C : H_x from H(x) = 0, defining the discontinuity surface 
% tspan : define the time span of the integration
% Initcond : the initial condition
% Modified by Peter Tang, 23rd/Feb/2022
function [tout,yout,yeout0,teout,yeout,ieout]=...
    Single_DS_Impacting_Hybrid_system_integration(A,R,C,InitCond,tspan,fs,equi_type)


% define integration timestep and timespan
tstart = tspan(1);
tfinal = tspan(2);
stepsize = 1/fs;
% define initial condition to start
y0 = InitCond;

% define the matrix of this linear differential system A
A = real(A);
% define the linear ode system
f_ls= @(t,y) A*y;
% if R(x) =x + w(x)*v(x) = (I + W)*x = R*x
%  while v(x) = C*A*x, w(x)_x =B, R = I + W = I + B*C*A; B = (R -
%  I)*A^(-1)*C';
B   =  (R -eye(size(R,1)))/A*C';

% acccoring to Nordmark's notations
A_s = (eye(size(A))- (B*C*A)/(C*A*B))*A;
% [V,D]=eig(A_s)
% diag(D)
%  ----determine the pseudo equilibrium x_pso on the boundary-------
% C * x_pso = 1; 
% C * A * x_pso = 0;
% A_s * x_pso= 0;

MM = [C;C*A;A_s];
b = [1;0;zeros(length(C),1)];
x_pso=MM\b;


% construct a transform matrix
% Trans = eye(length(C));
% Trans(C>0,:) = []; % get a (n-1) x n tranform matrix
% %
% M = A_s;
% % ROW_M  = C*M; % to check the vadality of the remaining equation
% COL_M = Trans*M(:,C>0);
% %
% M(C>0,:) =[];
% M(:,C>0) =[];
% Ans = -M\COL_M;
% % substitute the sub-system's solution to the remaining eqaution
% x_pso = (C'+ Trans'*Ans);
%------------------------------------END----------------------------------%
% define the sliding vector field
f_stck = @(t,y) A_s*(y-x_pso);

%-------------------------------BEGIN INTEGRATION-------------------------%

% flag is true when the flap stick to the boundary
BeStuck=false;

% start calculation
tout = tstart;
yout = y0.';
teout = [];
yeout = [];
yeout0 = [];
ieout = [];

% initialization for different contact mode
Event_flag=[];
% check if start at the boundary
if abs(C*y0-1)<1e-3 && C*A*y0>0
    [tstart,y0,BeStuck,tout,yout]=Post_impact(y0,A,R,C,tstart,tfinal,tout,yout,fs);
end
% define options
refine = 1;
options = odeset('Events',@(t,y) DetectingContact_events(t,y,C,equi_type),...
    'RelTol',1e-12,'AbsTol',1e-12,'Refine',refine);
opsl_stuck=odeset(options,'Events',@(t,y) ...
    unstickEvent(t,y,A,C));

while tout(end)<tfinal
    % define out put time serias
    try
        timespan=[tstart:1/fs:tfinal];
    catch ME
        keyboard
    end
    %     timespan=[tstart tfinal];
    
    % exit case I
    if length(timespan)<2
        % last time step
        break;
    end
    
    if BeStuck==false
        % Solve until the first terminal event: hit the boundary.
        [t,y,te,ye,ie] = ode45(@(t,y) f_ls(t,y),...
            timespan,y0,options);
        %         options = odeset(options,'InitialStep',t(end)-t(end-1));
        options = odeset(options,'InitialStep',1e-4*stepsize,'MaxStep',1e-2*stepsize);
        %         disp('hit the boundary one time')
        Event_flag='hit boundary';
        
        % In case for the small wrong initial stepsize
        check_discre=norm(y0-y(end,:));
        while check_discre<1e-6
            % In this case, failure may be due to static equilbrium point
            % in phase space(fixed point)
            t_flow= A*y0;
            if norm(t_flow)<1e-6
                disp('fixed point')
                BeStuck=true;
                break
            end
            
            NewInitialStep=2*options.InitialStep;
            sprintf('amplify initial stepsize')
            options = odeset(options,'InitialStep',NewInitialStep);
            [t,y,te,ye,ie] = ode45(@(t,y) f_ls(t,y)...
                ,timespan,y0,options);
            break;
        end
        %
        
    else
        % sliding contact
        opsl_stuck=odeset(opsl_stuck,'InitialStep',1e-2*stepsize,'MaxStep',stepsize);
        [t,y,te,ye,ie] = ode45(f_stck,timespan,y0,opsl_stuck);
        Event_flag='Stuck event over'
    end
    
    % Accumulate output.  This could be passed out as output arguments.
    nt    = length(t);
  
    % output history
    tout  = [tout; t(2:nt)];
    yout  = [yout; real(y(2:nt,:))];
    
    % output events
    teout = [teout; te];
    yeout = [yeout; ye];
    ieout = [ieout; ie];
    if size(y0,1)>1
        yeout0 = [yeout0; y0'];
    else
        yeout0 = [yeout0; y0];
    end
    %evaluate new initial condition for next integration step if possible
    if isempty(ie)
        % means that no event detected and the integration finished till
        % the time destination, so just end the intgration and return
        % result to the mainfunction
        break
    elseif BeStuck == true
        % if previous integration is sliding, then it shows the sliding
        % stoped, continue integration
        y0    = y(nt,:);
        tstart = t(nt);
        BeStuck = false; % stuck over
        continue
    end
 
    %------------------- else applay impact map--------------------------%
    t0=t(end);
    x=y(end,:)';
    [tstart,y0,BeStuck,tout,yout]=Post_impact(x,A,R,C,t0,tfinal,tout,yout,fs);

end





%--------------------------Detecting impact -----------------------------%
function [value,isterminal,direction] = DetectingContact_events(t,y,C,equi_type)
% Locate the time when flap hit through boundary in both directions and
% stop integration when isterminal ==1 Events: switch boundary; zero heave
% velocity; zero pitch velocity; zero flap velocity
value = C*y-equi_type*1;                %  detect switching point
isterminal = 1;                              %  stop the integration
direction  = -1;                             %  mark with event funct-
%  ion adding direction
%---------------------Detecting the taking off---------------------------%
function [value,isterminal,direction] = unstickEvent(t,y,A,C)
% detect the unstick event when the acceleration becomes bigger than 0
value =C*A*A*y;        %  detect event
isterminal = 1;                                   %  stop the integration
direction  = 1;                                   %  mark with event funct-
%  ion adding direction

%----------------- Post dynamics after contanct the \Sigma---------------%
function [tstart,y0,BeStuck,tout,yout]=Post_impact(x,A,R,C,t0,tfinal,tout,yout,fs)
Fx=A*x;
vx=C*Fx;
ax=C*A*A*x;
if vx>0 && t0<tfinal
    %  this case is abnormal and need further attention
    %         keyboard
    tstart=t0;
    y0=x;
    BeStuck=0;
    fprintf('positive vx %e at time %f !\n',vx,t0)
else % negative vx
    % end this chattering sequence when the velocity meet the threshhold
    Wx=(R -eye(size(R,1)))/A*C';
    rx = abs(1+C*A*Wx);
    if abs(vx)<1e-4 && ax<0 && t0<tfinal
        % apply chattering mapping
        % I: Nordmark's method to end the chattering
        xstar=x+(1/(1-rx) * ( 2 * Fx * rx/ax + Wx ))*vx;
        tstar =t0+1/(1-rx) * ( 2*rx/ax ) * vx;
        % II: My method to end the chattering
        % construct a transform matrix
%         Trans = eye(length(C));
%         Trans(C>0,:) = []; % get a (n-1) x n tranform matrix
%         %
%         n_v      = (C*A)';
%         n_v(C>0) = [];
%         xstar =x;
%         while abs(vx)> 1e-9
%         n_dis = abs(vx/norm(n_v));  
%         xstar = xstar - Trans'*n_dis*n_v*sign(vx);
%         vx=C*A*xstar;
%         end
      
        
        if abs(C*xstar-1)>0
            % this map occurrs when hitting the boundary but it's wrong
            % when the flag freedom penetrates the boundary, so we need fix
            % this mannuly with physical meaning
%             xstar(C>0)=1;
        end
        % if there is big jump in state variables, take some modification
        % for the map declare ending with the stick status
        BeStuck=true;
        
        % append the map resulted state to the output
        tout=[tout;tstar];
        yout=[yout;xstar'];
        
        % new initial condition for next integration
        tstart=tstar;
        y0=xstar;
        vx=C*A*y0;
    else % normal impact
        %Hit the boundary and rebounce with a new velocity with reverse
        %direction R(x) = (I + W)*x
        y0=R*x;
        tstart=t0+0.1/fs;
        tout=[tout;tstart];
        yout=[yout;y0'];
        BeStuck =false;
    end
end

