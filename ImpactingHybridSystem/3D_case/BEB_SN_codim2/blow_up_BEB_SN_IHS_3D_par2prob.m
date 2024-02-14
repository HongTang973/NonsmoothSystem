%< par to the problem 
%> prob contains the ode function: f_ode = prob.odes
function prob = blow_up_BEB_SN_IHS_3D_par2prob(par, equi_type)

%< use the function to form the matrixes 
[A,B,C,R]  = par2NForm_Lienard(par);

%> since the par(4) = mu will be in [0 mu_crit] but the equlibrium is
%> addmissible

%< the case mu = 0 should be futher investigated
% prob.equi_type          = -sign(par(4));
prob.equi_type          = equi_type;
%> get the initial condition as the fixed point from a 
Par_selected = ...
  [-0.1,0.2,-0.525,1.78192697901049+0.0955685096979647,1.6,6.674]';
[A1,B1,C1,R1]           = par2NForm_Lienard(Par_selected);
prob.IC                 = IC_generator(Par_selected(6),R1,A1,C1,prob.equi_type);
%>

prob.odes_Fcn  = @(t,y) A*y;

%> the period is around 6 second, so run the simulation for 100 cycles
prob.tspan  = [0 8000];
%>
prob.fs    = 20;

%> 
 prob.C     = C;
%> define  the function to get the accelaration in the observed coordinate
prob.ax_Fcn = @(t,y) C*A*A*y;
prob.Wx_Fcn = @(y) -B;
prob.rx_Fcn = @(x) par(4)-1;

prob.Impact_map  = @(x) x - B*C*A*x; 

prob.Post_impact_Fcn = @Shilnikov_homoclinic_chaos_post_impact;

%% define options
refine = 1;
prob.options = odeset('Events',@(t,y) ...
    Shilnikov_homoclinic_chaos_DetectingCross_events(t, y , prob),...
    'RelTol',1e-12,'AbsTol',1e-12,'Refine',refine);

prob.opsl_stuck=odeset(prob.options,'Events',@(t,y) ...
     Shilnikov_homoclinic_chaos_unstickEvent(t,y, prob));

end