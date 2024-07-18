%< par to the problem 
%> prob contains the ode function: f_ode = prob.odes
function prob = full_BEB_SN_IHS_3D_par2prob(par)

%< use the function to form the matrixes 
[A,B,C,R]  = par2NForm_DummyVar(par);
mu         = par(7);
M          = [0;0;-1];
% C*A^(-1)*M = -1/d;  
%> since the par(4) = mu will be in [0 mu_crit] but the equlibrium is
%> addmissible

%< The H(x)     = C*x = x(1) = 0;
prob.equi_type          = 0;
%> get the initial condition as the fixed point from a 
Par_selected = ...
  [-0.1,0.2,-0.525,1.78192697901049+0.0955685096979647,1.6,6.674]';
[A1,B1,C1,R1]      = par2NForm_Lienard(Par_selected);
IC                 = IC_generator(Par_selected(6),R1,A1,C1,1);
%>
IC(1) = 0;
prob.IC = sign(C*A*IC)*IC;
%> f(t,x) = A(mu,eta)x  + mu*M;
f_ls           = @(t,y) A*y + M*mu;
prob.odes_Fcn  = f_ls;

%> the period is around 6 second, so run the simulation for 100 cycles
prob.tspan  = [0 8000];
%>
prob.fs    = 20;

%> 
 prob.C     = C;
%> define  the function to get the accelaration in the observed coordinate
prob.ax_Fcn = @(t,y) C*f_ls(0,f_ls(0,y));
prob.Wx_Fcn = @(y) -B;
prob.rx_Fcn = @(x) par(4)+par(8)-1;
prob.bx_Fcn = @(y) C*A*B;
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