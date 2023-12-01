function [tout,yout,teout,yeout,yeout0,ieout,amplitude]=PWSC_int(A_1,InitCond,state0,tspan,fs,mu,flag)
%*************************#####************************* #
% checked and modified by Peter Tang /11th,3,2021       *#
%#************************#####************************* #
% gap:  backlash amount
% mu: scaling factor
% flag: return when there is event
% flag =0;
gap =0;
% define integration timestep and timespan
tstart = tspan(1);
tfinal = tspan(2);
% dt_scale = 1000/fs;
amplitude=[];
% define initial condition to start
y0 = InitCond;
A_matrix =A_1;
AA = A_1*A_1;
f_ = -A_1*state0;
%define options
refine = 1;
options = odeset('Events',@(t,y) DetectingZero_events(t,y),'RelTol',1e-9,'AbsTol',1e-9,...
    'Refine',refine);
opsl_stuck=odeset(options,'Events',@(t,y) ...
    unstickEvent(t,y,AA,state0));

% define the reset map
a11=A_1(1,1);a12=A_1(1,2);a13=A_1(1,3);
V_l3 = [0;a13;-a12]/sqrt(a12^2+a13^2);
P_l3 = [0;a12; a13]/sqrt(a12^2+a13^2);
r_ = @(y,a12,a13,f1) abs(a12*y(2)+a13*y(3)+f1)/sqrt(a12^2+a13^2);
z_ = @(y,V_l3) V_l3'*y;
%
% phi_1 = @(r_,z_) (1+0.9)*r_;
% phi_2 = @(r_,z_) 0*z_;
% %
% R = @(y,a12,a13,f1) y + phi_1(r_(y,a12,a13,f1),z_(y,V_l3))*P_l3 + phi_2(r_(y,a12,a13,f1),z_(y,V_l3))*V_l3;
% function during region I
f_ls= @(t,y) A_matrix*(y-state0);
% functions needed for chattering mapping and the stuck case
pr=0.8; %coef of restitution
pz = 4.6;
ax_= @(t,y) [1 0 0]*AA*(y-state0);
vx_= @(t,y)  [1 0 0]*A_matrix*(y-state0); 
N = sqrt(a12^2+a13^2);
% MW = @(y)  y +(-(1+pr)*P_l3*A_1(1,:)/N+pz*V_l3*A_1(1,:)/N)*(y-state0);
Wx_=@(t,y) (-(1+pr)*[0 a12 a13]'+pz*[0 a13 -a12]')./(a12^2+a13^2);
z_left = @(y) -(f_(1)+a12*y)/a13;
% alternative function to use when stuck in sliding contact
f_stck = @(t,y) f_ls(t,y) + (ax_(t,y)/(1+pr))*Wx_(t,y);
%------------------------------------------------------------------------%
% the stability of the sticking set
A_s = (eye(3) - (Wx_(0,0)*[1 0 0]*A_1)/([1 0 0]*A_1*Wx_(0,0)))*A_1;
[Vs,Ds] = eig(A_s);
% pseudo equilibrium
% zp = (A_s(1,1)*A_s(2,2)-A_s(2,1)*A_s(1,2))/(A_s(2,3)*A_s(1,2)-A_s(2,2)*A_s(1,3))
% yp = z_left(zp)
% yp = (A_s(1,1)*A_s(2,3)-A_s(1,3)*A_s(2,1))/(A_s(1,3)*A_s(2,2)-A_s(2,3)*A_s(1,2));
% zp = z_left(yp);
% flag is true when the flap get steady contact with the constrain boundary
BeStuck=false;


% start calculation
tout = tstart;
yout = y0.';
teout = [];
yeout = [];
yeout0= [];
ieout = [];

% check if start at the boundary
if abs(y0(1))<1e-3 && ax_(0,y0)<0
    [tstart,y0,BeStuck,tout,yout]=Post_impact(y0,f_ls,pr,vx_,ax_,Wx_,z_left,tstart,tfinal,tout,yout,@impact_map,fs,gap,mu);
end

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
    if BeStuck==false
        % Solve until the first terminal event: hit the boundary.
        [t,y,te,ye,ie] = ode45(@(t,y) f_ls(t,y),...
            timespan,y0,options);
        %         options = odeset(options,'InitialStep',t(end)-t(end-1));
        options = odeset(options,'InitialStep',1e-4,'MaxStep',1e-3);
        %         disp('hit the boundary one time')
        Event_flag='hit boundary';
        
        % In case for the small wrong initial stepsize
        check_discre=norm(y0-y(end,:));
        while check_discre<1e-6
            % In this case, failure may be due to static equilbrium point
            % in phase space(fixed point)
            t_flow= A_matrix*y0+state0(:,1);
            if norm(t_flow)<1e-6
                disp('fixed point')
                BeStuck=true;
                break
            end
            
            NewInitialStep=2*options.InitialStep;
            sprintf('amplify initial stepsize')
            %             keyboard
            options = odeset(options,'InitialStep',NewInitialStep);
            [t,y,te,ye,ie] = ode45(@(t,y) f_ls(t,y)...
                ,timespan,y0,options);
            break;
        end
        %
        
    else
        % sliding contact
        opsl_stuck=odeset(opsl_stuck,'InitialStep',1e-4,'MaxStep',1e-2);
        f_stck(tstart,y0);
        [t,y,te,ye,ie] = ode45(f_stck,timespan,y0,opsl_stuck);
        f_stck(t,y(end,:)');
        Event_flag='Stuck event over'
    end
    
    
    r1 = r_(y0,a12,a13,f_(1));
    try 
    z1 = z_(y0,V_l3);
    catch ME
        keyboard
    end
    if ~isempty(ye)
        r2 = r_(ye,a12,a13,f_(1));
        z2 = z_(ye',V_l3);
        amplitude=[amplitude,[r2/r1;z2/z1]];
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
    yeout0 = [yeout0; y0'];
    yeout = [yeout; ye];
    ieout = [ieout; ie];
    if flag % when there is need to halt at first event
        break
    end
    %     if size(y0,1)>1
    %     yeout0 = [yeout0; y0'];
    %     else
    %         yeout0 = [yeout0; y0];
    %     end
    %evaluate new initial condition for next integration step if possible
    if isempty(ie)
        % means that no event detected and the integration finished till
        % the time destination, so just end the intgration and return
        % result to the mainfunction
%         vx=vx_(t0,y(end,:)');
%         vx=vx_(t0,[0 yp zp]');
%         ax=ax_(t0,[0 yp zp]');
        
        break
    elseif BeStuck == true
        % if previous integration is sliding, then it shows the sliding
        % stoped, continue integration
        y0    = y(nt,:)';
        tstart = t(nt);
        BeStuck = false; % stuck over
        continue
    end
    %    dsiP(impact_map)
    %------------------- else applay impact map--------------------------%
    t0=t(end);
    x=y(end,:)';
    
    [tstart,y0,BeStuck,tout,yout]=Post_impact(x,f_ls,pr,vx_,ax_,Wx_,z_left,t0,tfinal,tout,yout,@impact_map,fs,gap,mu);
    %     A_1*ye'
    % Set the new initial conditions for next integration
    
    %
 
end





%------------------------------------------------------------------------%
%------------------------------------------------------------------------%


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

%
function [value,isterminal,direction] = unstickEvent(t,y,AA,state0)
% detect the unstick event when the acceleration becomes bigger than 0
Fx=AA*(y-state0);
 value =([1 0 0]*Fx);        %  detect event
isterminal = 1;                                   %  stop the integration
direction  = 1;


% function reset map
function y0=impact_map(vx,Wx,x)
y0 = x + Wx*vx;


function [tstart,y0,BeStuck,tout,yout]=Post_impact(x,f_ls,rx,vx_,ax_,Wx_,z_left,t0,tfinal,tout,yout,impact_map,fs,gap,mu)
vx=vx_(t0,x);
Fx=f_ls(t0,x);
ax=ax_(t0,x);
Wx=Wx_(t0,x);
if vx>0 && t0<tfinal
    %  this case is abnormal and need further attention
    %         keyboard
    tstart=t0;
    y0=x;
    BeStuck=0;
    fprintf('positive vx %e at time %f !\n',vx,t0)
else % negative vx
    
    
    % end this chattering sequence when the velocity meet the threshhold
    if abs(vx)<mu*1e-5 && ax<0 && t0<tfinal
        % apply chattering mapping
        
        if rx==1
            xstar=x+( 1/(1-rx+exp(-0.5/fs)) ) * ( 2 * Fx * rx/ax + Wx ) *vx;
            tstar =t0+(1/(1-rx+exp(-0.5/fs)) ) * ( 2*rx/ax ) * vx;
        else
            xstar=x+(1/(1-rx) * ( 2 * Fx * rx/ax + Wx ))*vx;
            tstar =t0+1/(1-rx) * ( 2*rx/ax ) * vx;
            
        end
        %         fprintf('chattering detected at vx %e time %f !\n ',vx,t0)
        %
        if abs(xstar(1))>gap
            % this map occurrs when hitting the boundary but it's wrong
            % when the flag freedom penetrates the boundary, so we need fix
            % this mannuly with physical meaning
            
            %             disp('wrong chattering map')
            xstar(1)=gap*sign(xstar(1));
            xstar(3) = z_left(xstar(2));
%             vx=vx_(t0,xstar)
        end
        % if there is big jump in state variables, take some modification
        % for the map
        %         data_h=yout(:,1);
        %         IndMin=find(diff(sign(diff(data_h)))>0)+1;   %获得局部最小值的位置
        %         IndMax=find(diff(sign(diff(data_h)))<0)+1;   %获得局部最大值的位置
        %         x_h1=[yout(IndMin(end),1),yout(IndMax(end),1)];
        %         data_alpha=yout(:,2);
        %         IndMin=find(diff(sign(diff(data_alpha)))>0)+1;   %获得局部最小值的位置
        %         IndMax=find(diff(sign(diff(data_alpha)))<0)+1;   %获得局部最大值的位置
        %         x_alpha1=[yout(IndMin(end),2),yout(IndMax(end),2)];
        %
        %         if xstar(1)>x_h1(1)&&xstar(1)<x_h1(2)&&xstar(2)>x_alpha1(1)&&xstar(2)<x_alpha1(2)
        %         else
        %             xstar(1:2)=0.5*[sum(x_h1),sum(x_alpha1)];
        %             disp('map modification II')
        %         end
        %         fprintf('Chattering complete at time %f !\n' , tstar );
        % declare ending with the stick status
        BeStuck=true;
        
        % append the map resulted state to the output
        tout=[tout;tstar];
        yout=[yout;xstar'];
        
        % new initial condition for next integration
        tstart=tstar;
        y0=xstar;
        
    else % normal impact
        %Hit the boundary and rebounce with a new velocity with reverse
        %direction
%         vx=vx_(t0,y0)
        y0=impact_map(vx,Wx,x);
%         y0=mu(x);
%         y0=impact_map(x-state0,mu,x);
%         vx=vx_(t0,y0)
        tstart=t0+0.1/fs;
        tout=[tout;tstart];
        yout=[yout;y0'];
        % delta_v=MHH(1:2,1:2)\[MHH(1,3);MHH(2,3)]*(1+res)*y0(6);
        % y0(6)=-res.*y0(6);
        % y0(4)=y0(4)+delta_v(1);
        % y0(5)=y0(5)+delta_v(2);
        %         disp('bounce one time')
        
        BeStuck =false;
    end
end


