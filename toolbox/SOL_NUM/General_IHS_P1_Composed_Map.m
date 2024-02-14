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
C       = prob.OB_C;
if isfield(prob, 'IPC_num_allowed')
    IPC_num_allowed = prob.IPC_num_allowed;
else
    IPC_num_allowed =1;
end
%> I ---- operator one: the flow evolution
timespan = 1.2*[0 T];

%> define the event function
EventFun = prob.efunc;
%>
refine  = 1;
%     EventFun(0,y0)
options = odeset('Events',@(t,y) EventFun(t,y),...
    'RelTol',1e-12,'AbsTol',1e-12,'Refine',refine);

%> this is a general differential equation which is determined by the
%> parameters of the system, and could be nonlinear 
flow_operator = prob.par_2_flow_operator(sys_par);

ResetMap      = prob.par_2_ResetMap(sys_par);


%% --------------- check the fixed point condition ---------------- 
if strcmp(prob.evol_seq,'phi_map')
    %> initial step: Solve until the first terminal event: hit the boundary.
    [t1,y1,te_0,ye,~] = ode45(@(t,y) flow_operator(t,y), timespan, x_p, options);
    %> the returning point
    x_00    = ResetMap( y1(end,:)' );
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
        [t1,y1,te,ye,ie] = ode45(@(t,y) flow_operator(t,y), timespan, IC, options);
        %> the returning point
        if ~isempty(ie)
            IC    = ResetMap( ye' );
            % Get the new run of flow  with new initial conditions
        else
            figure; plot(t1,y1(:,C>0));
            error('No impact happens!')
        end
    end
else %>impose reset map first
    for i = 1:IPC_num_allowed
        %> initial step: Solve until the first terminal event: hit the boundary.
        [t1,y1,te,ye,ie] = ode45(@(t,y) flow_operator(t,y), timespan, ResetMap(IC), options);
        %> the returning point
        if ~isempty(ie)
            IC    = ye';
            % Get the new run of flow  with new initial conditions
        else
            figure; plot(t1,y1(:,C>0));
            error('No impact happens!')
        end
    end
end
%------------------- return the image of the map after IPC_num_allowed
%   times of evolution   ---------------------------------------------
Maped_state = [IC;te];








