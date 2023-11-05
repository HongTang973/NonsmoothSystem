
% function [u,iter] = codim1_PC(ZeroFunctions,par,ind_p,ind_x,p_span,x_span,dir)

% Parameter Analysis of the found period one orbit

%> functionality added on 16/10/2023  -- Peter
%> --- the ability to detect the PD or SN along the continuation curve
%> every time when a new par is created, a tag of classification should be
%> given



% Self coded by knowledg from the book
% Introduction to Numerical Continuation Methods --- Allgower

% ind_p:            specify the varying parameter
% ind_x:            specify the solution state
% ZeroFunctions:    function handle of zero funcitons of OP
% par:              give the para meter vector
% p_span:           specify the range of parameter to be continued
% dir:              specify the direction of the solution branch
function varargout = codim1_PC(prob)

% function [u,iter] = codim1_PC(prob)
%
ind_p               = prob.ind_p;
ind_x               = prob.ind_x;
ZeroFunctions       = prob.ZeroFunctions;
par                 = prob.par;
% keys                = prob.keys;
p_span              = prob.p_span;
x_span              = prob.x_span;
dir                 = prob.dir;
%> added on 26/10/2023
bifur_detct         = prob.bifur_detct;
NumOfEvents         = prob.Events.NumOfEvents;
MonitorFunctions    = prob.Events.MonitorFunctions;
Events_tag          = prob.Events.tag;
lab =[];
lab                 = [lab;'SP']; %> 0 for plain points; 1 for SN points; -1 for the PD points
if bifur_detct
    Monitor_fval    = [];
    Monitor_fval        =  [ Monitor_fval, MonitorFunctions(prob,par)];
end

%
if isvector(par)
    par = par;
elseif isstruct(par)
    par = cell2mat(struct2cell(par));
else
    warndlg('Wrong input in the codim1_PC function!')
end
% now par must be a vector now
% arc length stepsize

if isfield(prob,'uplim_step')
    uplim_step = prob.uplim_step;
else
    uplim_step = 100;
end

if isfield(prob,'co_dh')
    default_h = prob.co_dh;
else
    default_h = 0.001;
end

dh = min(default_h,(p_span(2)-p_span(1))/200);

% DETERMINE THE LOCATION of the parameter's partial derivative
index = [ind_p,ind_x];

if  size(x_span,1)~= length(ind_x) || size(p_span,1)~= length(ind_p)
    warndlg('Wrong length match of p and x.')
end

% searching direction of parameter space
% dir = dir;


%   **************** Codimension one ******************************    %


% step 1: set the closed form zero equations
% H(x,a) = 0
% u_0 = [x,a];
% H(u_0) = 0;

% step 2:  using the local manifold information to span the solution branch
%  Jacobian matrix : H'(u_0)
%  Rank(H'(u_0)) = n;  the one-dimensional kernel ker(H')
%  if u(a) denotes the solution branch
%               H'(u) u'(a) = 0
%  use the condition of det([H'; u']) > 0 to choose the right tangent
%  vector


% step 3: stepsize of prediction by the arclength
% ds =  norm(du_da) da



% step 4: Euler predictor


% step 5: corrector step

%  ---------------------       ----------------------------      -------- %

% ##### Step 1: define zero functions

% scalar function 1D searching



% check the par is a vector
if size(par,2)==1;  else;  par = par.';  end

% initial condition of the continuation system

u = par; % container of the  parameter track

iter = [];



[prob_det] = ZeroFunctions (par);

if norm(prob_det) > 5e-9
    fprintf( 'Not zero parameter point! \n' )
    % correct the initial condition
    [u,~,~,~] = corrector(ZeroFunctions,par, par, ind_p, ind_x,0,0);
    par                 =  u(:,end);
    if bifur_detct
    Monitor_fval        =  [ Monitor_fval, MonitorFunctions(prob,u(:,end))];
    lab                 =  [lab;'SP'];
    end
    fprintf('The SP is relocated at curve with evalf f = %g \n',ZeroFunctions (par))
end
%

% ##### Step 2: get the Jacobian
% par = {'U', 'r', 'damp_1','damp_2', 'damp_3', 'T'}



[t_,~] = get_jacobian(ZeroFunctions, par, index);
% get the tangent vector induced by the Jacobian
% t_ = null(Jacob);
%
% if sign(det([Jacob;t_']))==1
%
% elseif sign(det([Jacob;t_']))==-1
%     t_ = -t_;
% else
%     t_ = 0;
% end

%>
trans_flag =  norm(t_);

reg = u(ind_p,end)>=p_span(1) && u(ind_p,end)<=p_span(2)...
    && all(u(ind_x,end)-x_span(:,1)>=0) && all(x_span(:,2)-u(ind_x,end)>=0);

new_par = u(:,end); dis2init = [1,1];
nstep   = 1;
while trans_flag && reg && nstep < uplim_step   
    % first predictor
    % using the tangent vector produced by the Jacobian
    new_par = Predictor(new_par,index,ind_p,t_,dh,dir,dis2init);
    % corrector: now the newly predicted point in parameter space, but this 
    % should be corrected by iteration
    num = 1;
    [u,iter,dis2init,iscirc,discre_err] = corrector(ZeroFunctions,u, new_par, ind_p, ind_x,num,iter);
    fprintf('Convergence to %3f with iter %d times with err %g! \n', u(ind_p,end),iter(end),discre_err);
    %>  added on 16/10/2023 ----- Event monitor 
    if bifur_detct
        %> get the monitor function values to check if the events are
        %> happening
        [Events_favls, S_US]  = MonitorFunctions(prob,u(:,end));
        %> check the sign change of the function
        if size(Monitor_fval,2) == 0
            bull_events_check = ones(size(Events_favls));
            %>
        else
            bull_events_check = sign(Events_favls).*sign(Monitor_fval(:,end));
        end
        %> give tag to the points
        if any(bull_events_check == 0) %> just record the event tag and continues to do the continuation
            %> this rarely happens  and nearly impossible when two zeros
            %> happen simutaneously
            event_tag     = Events_tag(bull_events_check == 0);
            if sum(bull_events_check == 0)>1;  error('Two events happening at the same time! \n'); end
            %> append the monitor function values in the buffing
            Monitor_fval  = [Monitor_fval, Events_favls];
            lab           = [lab; event_tag];
            
        elseif any(bull_events_check < 0) %> puase for the root finding
            %> get the last step resuts as tmp
            %> u,iter,dis2init,iscirc,discre_err
            u_last_step          = u(:,end);                   u(:,end) =[];
            iter_last_step       = iter(end);                 iter(end) =[];
            dis2init_last_step   = dis2init(:,end);     dis2init(:,end) =[];
            %             iscirc_last_step     = iscirc(end);             iscirc(end) =[];
            discre_err_last_step = discre_err(end);     discre_err(end) =[];
            %> to locate the bifurcation point
            event_tag = Events_tag(bull_events_check < 0);
            if sum(bull_events_check < 0)>1;  keyboard;error('Two events happening at the same time! \n'); end
            %> the intermediate value should be find to locate the exact
            %> bifur point
            ind = bull_events_check < 0; %> this is the unit vector to select the function
            Event_fun = @(par) (ind'* MonitorFunctions(prob, par))^2;
            Bifur_zero_Fcns = @(par) [ZeroFunctions(par);Event_fun(par)];
            %> use the newton method in the corrector step to sove the zero
            %> solution
            fprintf('The event %s is detected! \n',event_tag{1});
            %> use the iteration process to get the critical point of PD/SN
            %------------------------------------------------------ %
            new_par     =  u_last_step;
            [prob_det]  =  Bifur_zero_Fcns (new_par);
            dis         =  norm(prob_det);
            num         = 1;
            
            while dis  > 1e-4
                [prob_det] =  Bifur_zero_Fcns (new_par);
                dis = norm(prob_det);
                [~,J] = get_jacobian( Bifur_zero_Fcns, new_par, index);
                alt = -real(J'*inv(J*J')*prob_det);
                num = num +1;
                if norm(alt)>1e1
                    keyboard
                    alt = alt./norm(alt)/100;
                end
                new_par(index) = new_par(index) + alt;  
            end
            %------------------------------------------------------ %
            %> record 
            discre_err = [discre_err, dis];
            iter       = [iter, num];
            u          = [u, new_par];
            dis2init   = [dis2init, dis2init_last_step];
            %             [u,iter,dis2init,iscirc,discre_err] = corrector(Bifur_zero_Fcns,u, u_last_step, ind_p, ind_x,1,iter);
            %> record this extra interpolated point with tag
            Monitor_fval  = [Monitor_fval, MonitorFunctions(prob,u(:,end))];
            lab           = [lab; event_tag];

            %> record the ending point before the interpolation
            u             = [u,        u_last_step];
            iter          = [iter,     iter_last_step];
            dis2init      = [dis2init, dis2init_last_step];
            %             iscirc        = [iscirc,   iscirc_last_step    ];
            discre_err    = [discre_err,discre_err_last_step];
            %>
            Monitor_fval  = [Monitor_fval, Events_favls];
            lab           = [lab; {S_US}];
            
        else %> no pause and continues to continuation
            event_tag = {S_US};
            %> append the monitor function values in the buffing
            Monitor_fval  = [Monitor_fval, Events_favls];
            lab           = [lab; event_tag];
        end
                   
    end
    %>
    [t_,~] = get_jacobian(ZeroFunctions, u(:,end), index);
    %>
    trans_flag =  norm(t_);
    %>
    reg = iscirc && u(ind_p,end)>=p_span(1) && u(ind_p,end)<=p_span(2)...
        && all(u(ind_x,end)-x_span(:,1)>=0) && all(x_span(:,2)-u(ind_x,end)>=0);
    %>
    nstep  = nstep + 1;
end

%> output:
if nargout >= 0; varargout{1} =  u;             end
if nargout >= 2; varargout{2}  = iter;          end
if nargout >= 3; varargout{3}  = lab;           end


%> ------------- built in funcitons -------------------------- %

function new_par = Predictor(new_par,index,ind_p,t_,dh,dir,dis2init)

if abs(t_(1))<1e-2  % turning point
    % the index =1, since we put the parameter in the first place when
    % we get our jacobian
    
    dh = max( 1e-3, min(0.02*dh, 0.1*abs(dis2init(1))));
    K = null(t_');
    co = [1;rand(length(t_)-1,1)];
    %> strategy 1
    co = co/norm(co);
    % new searching direction
    t_ = [t_,K]*co;
    
elseif  abs(dis2init(1))<2*dh && sign(dis2init(1)*t_(1))<0 && dis2init(2)<1e-3
    % the second condition is the loop returning conditon to prevent the
    % low starting speed at first steps
    % dh = min(dh*0.1,0.5*abs(dis2init(1))*new_par(ind_p));
    dh = min(dh*0.1, 0.001);
    
else
    % plain predictor step
    
end
new_par(index) = new_par(index) + dir * dh* t_;


%%
function [u,iter,dis2init,iscirc,discre_err] = corrector(ZeroFunctions,u, new_par, ind_p, ind_x,num,iter)
index   = [ind_p,ind_x];
par     = u(:,1); % the first point
while 1
    [prob_det] = ZeroFunctions (new_par);

    [~,J] = get_jacobian(ZeroFunctions, new_par, index);

    alt = -real(J'*inv(J*J')*prob_det);
    
    new_par(index) = new_par(index) + alt;

    num = num + 1;


    % check convergence
    discre_err = norm(prob_det);
    if  discre_err <= 1e-12 || num >=100

        u = [u,new_par];
        iter = [iter,num];


        % incase a closed circle solution
        dis2init = [(new_par(ind_p)-par(ind_p))/(1+norm(par(ind_p)));
            norm(new_par(ind_x)-par(ind_x))/norm(par(ind_x))] ;

        iscirc = norm(new_par(index)-par(index))/(1+norm(par(index)))>1e-12;

        if num >=100
            fprintf('Iteration times limited!! with %d',norm(prob_det))
        end

        %
        break;
    end
end

