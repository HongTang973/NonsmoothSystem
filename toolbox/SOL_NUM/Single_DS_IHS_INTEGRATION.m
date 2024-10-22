%  Numerical integration of the Impacting Hybrid system with:
% the A matrix: F_ls = A*y  defing the flow of the linear differential
% equation 
% R : the reset map matrix, this means the flow is continuous 
% C : H_x from H(x) = 0, defining the discontinuity surface 
% tspan : define the time span of the integration
% Initcond : the initial condition
% Modified by Peter Tang, 23rd/Feb/2022
function [tout,yout,yeout0,teout,yeout,ieout,BeStuck]=...
    Single_DS_IHS_INTEGRATION(prob)
%> input: prob

%> prob contains the ode function: f_ode = prob.odes
f_odes      = prob.odes_Fcn;
%> 
InitCond    = prob.IC;
%>
tspan       = prob.tspan;
%>
fs          = prob.fs;
%> 
C           = prob.C;
%> define the acceleration 
% functions needed for chattering mapping and the stuck case
%coef of restitution rx



%> ax_= @(t,y,p) C*f_odes_jacobian(y)*f_odes(t, y, p);
ax_         =   prob.ax_Fcn;

%>  Wx_=@(t,y) -(1+rx)*[0 1 0 0 0]';
Wx_         =   prob.Wx_Fcn;

rx_         =   prob.rx_Fcn;

bx_         =   prob.bx_Fcn;
% alternative function to use when stuck in sliding contact 
f_stck      = @(t,y) f_odes(t,y)  +  ( ax_(t,y)/ bx_(y) ) .* Wx_(y);

% f_stck          =   prob.F_stick_Fcn;

%> define the impact map
Post_impact     = prob.Post_impact_Fcn;

% define integration timestep and timespan
tstart = tspan(1);
tfinal = tspan(2);
stepsize = 1/fs;
% define initial condition to start
y0 = InitCond;

%-------------------------------BEGIN INTEGRATION-------------------------%

% flag is true when the flap stick to the boundary
BeStuck=false;

% start calculation
tout    = tstart;
yout    = y0.';
teout   = [];
yeout   = [];
yeout0  = [];
ieout   = [];

% initialization for different contact mode
Event_flag=[];

% check if start at the boundary
if abs(C*y0-prob.equi_type)<1e-6 && C*f_odes(0,y0)<=1e-9
    [tstart,y0,BeStuck,tout,yout]=Post_impact(prob, y0,tout,yout);
end
% post_impact(prob,x,tout,yout) -- syntac
% define options

options         = prob.options;

opsl_stuck      = prob.opsl_stuck;

while tout(end)<tfinal
    % define out put time serias
    try
        timespan=[tstart:1/fs:tfinal];
%           timespan=[tstart tfinal];
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
        [t,y,te,ye,ie] = ode45(@(t,y) f_odes(t,y),timespan,y0,options);
%         options = odeset(options,'InitialStep',1e-2*stepsize,'MaxStep',stepsize);
        %         disp('hit the boundary one time')

        Event_flag='x0';
        
        % In case for the small wrong initial stepsize
        check_discre=norm(y0-y(end,:));
        while check_discre<1e-6
            % In this case, failure may be due to static equilbrium point
            % in phase space(fixed point)
            try
            t_flow= f_odes(0,y0);
            catch ME
                keyboard
            end
            if norm(t_flow)<1e-6
                % disp('fixed point')
                BeStuck=true;
                break
            end
            
            NewInitialStep=2*options.InitialStep;
            sprintf('amplify initial stepsize')
            options = odeset(options,'InitialStep',NewInitialStep);
            [t,y,te,ye,ie] = ode45(@(t,y) f_odes(t,y)...
                ,timespan,y0,options);
            break;
        end
        %
        
    else
        % sliding contact
        opsl_stuck=odeset(opsl_stuck,'InitialStep',1e-2*stepsize,'MaxStep',1000*stepsize);
        [t,y,te,ye,ie] = ode45(@(t,y) f_stck(t,y),timespan,y0,opsl_stuck);
        Event_flag='Stuck event over';
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
    if isempty(ie) || abs(tout(end) - timespan(end))<10*eps
        % means that no event detected and the integration finished till
        % the time destination, so just end the intgration and return
        % result to the mainfunction
        break
    elseif BeStuck == true
        % if previous integration is sliding, then it shows the sliding
        % stoped, continue integration
        y0    = y(nt,:)';
        tstart = t(nt);
        if ie == 2; break; end
        BeStuck = false; % stuck over
        continue
    end
 
    %------------------- else applay impact map--------------------------%
    t0=t(end);
    x=y(end,:)';
    [tstart,y0,BeStuck,tout,yout]=Post_impact(prob,x,tout,yout);

end








