function Maped_state = General_IHS_P1_Composed_Map(prob, sys_vec)
%> double mode: determined by prob.evol_seq
%> mode 1: map first -> free flight
%> mode 2: free flight first -> map ensues

%> single impact composed map evaluation
%< input:
%> IC：- fixed point with perturbation: x_00
%> IPC_num_allowed : Impact times -- default as 1
%> R : the reset map
%> T : the guessed upper bond of the flight time
%> C : the observation vector, normally being a unit vector with one unit
%> element

% sys_vec = [x_p; T_p; sys_par; IC; T];
x_p     = sys_vec(prob.sys_vec_index.xp_index);
T_p     = sys_vec(prob.sys_vec_index.Tp_index);
sys_par = sys_vec(prob.sys_vec_index.sys_par_index);
IC      = sys_vec(prob.sys_vec_index.IC_index);
T       = sys_vec(prob.sys_vec_index.T_simu_index);
C       = prob.C;
if isfield(prob, 'IPC_num_allowed')
    IPC_num_allowed = prob.IPC_num_allowed;
else
    IPC_num_allowed =1;
end
%> I ---- operator one: the flow evolution
timespan = 1.2*[0 T];


%> this is a general differential equation which is determined by the
%> parameters of the system, and could be nonlinear
flow_operator = prob.par_2_flow_operator(sys_par);
prob.odes_Fcn = flow_operator;
ResetMap      = prob.par_2_ResetMap(sys_par);

%> define the event function --- modified by Peter 1/8/2024 to consider
%non-invariant monitoring function
EventFun = @(t,y) prob.efunc(t,y,prob);
%>
refine  = 1;
%     EventFun(0,y0)
options = odeset('Events',@(t,y) EventFun(t,y),...
    'RelTol',1e-12,'AbsTol',1e-12,'Refine',refine);

%% --------------- check the fixed point condition ----------------
if strcmp(prob.evol_seq,'phi_map')
    %> initial step: Solve until the first terminal event: hit the boundary.
    [t1,y1,te_0,ye,~] = ode45(@(t,y) flow_operator(t,y), timespan, x_p, options);
    %> the returning point
    x_00    = ResetMap( y1(end,:)' );
elseif strcmp(prob.evol_seq,'phi_map_phi') %> poincare section at v =0
    %> initial step: Solve until the first terminal event: hit the boundary.
    [t1,y1,te_0,ye,~] = ode45(@(t,y) flow_operator(t,y), timespan, x_p, options);
    %> the returning point
    x_00    = ResetMap( y1(end,:)' );
    [t1,y1,te_0,ye,ie] = ode45(@(t,y) flow_operator(t,y), timespan, x_00, options);
    %> the returning point at section v =0
    x_00    = ye(ie==2,:)';
else
    %> initial step: Solve until the first terminal event: hit the boundary.
    [t1,y1,te_0,ye,~] = ode45(@(t,y) flow_operator(t,y), timespan, ResetMap(x_p), options);
    %> the returning point
    x_00    = y1(end,:)';
end

%> check if this already formulates the circle
if strcmp(prob.call_type,'jacobian')
    if norm(x_00 - x_p)/norm(x_p) > 5e-6
        figure; plot(t1,y1(:,C>0));
        fprintf('The x_p is not the starting point of a LCO! \n')
    end
else
    %
end
% -----------------------------------------------------------------

%> evolute the system with allowed impact times
if strcmp(prob.evol_seq,'phi_map')
    for i = 1:IPC_num_allowed
        %> initial step: Solve until the first terminal event: hit the boundary.
        [t1,y1,te_0,ye,ie] = ode45(@(t,y) flow_operator(t,y), timespan, IC, options);
        %> the returning point
        if ~isempty(ie)
            IC    = ResetMap( y1(end,:)' );
            % Get the new run of flow  with new initial conditions
        else
            figure; plot(t1,y1(:,C>0));
            error('No impact happens!')
        end
    end
    te = te_0(ie==1);
elseif strcmp(prob.evol_seq,'phi_map_phi') %> poincare section at v =0
    for i = 1:IPC_num_allowed
        %> initial step: Solve until the first terminal event: hit the boundary.
        [t1,y1,te_0,ye,~] = ode45(@(t,y) flow_operator(t,y), timespan, IC, options);
        %> the returning point
        x_00    = ResetMap( y1(end,:)' );
        [t1,y1,te_1,ye,ie] = ode45(@(t,y) flow_operator(t,y), timespan, x_00, options);
        %> the returning point at section v =0
        IC    = ye(ie==2,:)';
    end
    te = te_0 + te_1(ie==2);
else %>impose reset map first
    for i = 1:IPC_num_allowed
        %> initial step: Solve until the first terminal event: hit the boundary.
        [t1,y1,te_0,ye,ie] = ode45(@(t,y) flow_operator(t,y), timespan, ResetMap(IC), options);
        %> the returning point
        if ~isempty(ie)
            IC    = y1(end,:)';
            % Get the new run of flow  with new initial conditions
        else
            figure; plot(t1,y1(:,C>0));
            error('No impact happens!')
        end
    end
    te = te_0(ie==1);
end
%------------------- return the image of the map after IPC_num_allowed
%   times of evolution   ---------------------------------------------
Maped_state = [IC;te];








