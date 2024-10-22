%< par to the problem 
%> prob contains the ode function: f_ode = prob.odes
function prob = full_BEB_SN_IHS_3D_par2prob_NL(par)

%< use the function to form the matrixes 
[A,B,C,R]  = par2NForm_DummyVar(par);
% A = [ a1 1 0;
%       a2 0 1; 
%       a3 0 0];

% B = [0; b2; b3];

% C = [1,0,0];
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
IC(1)   = 0*abs(mu);
prob.IC = sign(C*A*IC)*IC;

% dyft_1 = ... + a1*mu + mu^2*(b11*y(1) +b12*y(2) + b13*y(3) ) + mu*( c11*y(1)^2 + c12*y(1)*y(2) +
% c13*y(1)*y(3) + c23*y(2)*y(3) + c22*y(2)^2 + c33*y(3)^2);
% dyft_2 = ... + a1*mu + mu^2*(b21*y(1) +b22*y(2) + b23*y(3) ) + mu*( d11*y(1)^2 + d12*y(1)*y(2) +
% d13*y(1)*y(3) + d23*y(2)*y(3) + d22*y(2)^2 + d33*y(3)^2);
% dyft_2 = ... + a1*mu + mu^2*(b31*y(1) +b32*y(2) + b33*y(3) ) + mu*( e11*y(1)^2 + e12*y(1)*y(2) +
% e13*y(1)*y(3) + e23*y(2)*y(3) + e22*y(2)^2 + e33*y(3)^2);

% in matrix form 
% NL_term = B00*[x1*x2; x2*x3; x3*x1] + C00*[x1^2; x2^2; x3^2]
% --- index = [2 3 1]
% -- bilinear term:     y.*y(index)
% -- quadratic term:    y.*y
index = [2 3 1];
%% this is the orginally proposed case 
scale_1     = 1;%-0.1; 
scale_2     = 0.3;% 0.01; 
scale_3     = 0.3;% 0.01; 
% from the formulae, we conclude that the c22 c23 c33 should be zero to get
% uniformly linear impact restitution law

A00 = scale_1* [0.25    0.1     0.2;
                0       -0.2    0.15;
                0.15    0.2     -0.1];

% B00*[x1*x2; x2*x3; x3*x1]
B00 = scale_2* [-0.6438    0   -0.5821;
                -0.2807    -0.3283  0.8103;
                -0.8866   -0.6487   0.3508];
% C00*[x1^2; x2^2; x3^2]
C00 = scale_3* [0.25*4     0       0;
                -0.5    -0.3   0.35;
                0.45   -0.6   0.1];

%% The proposed case is choosen after discussion with David
scale_1     = 1;%-0.1; 
scale_2     = 1;% 0.01; 
scale_3     = 1;% 0.01; 
% from the formulae, we conclude that the c22 c23 c33 should be zero to get
% uniformly linear impact restitution law

A00 = scale_1* [0.5    0    0;
                0    0    0;
                0    0    0];

% B00*[x1*x2; x2*x3; x3*x1]
B00 = scale_2* [-1    0   -1;
                0    0   0;
                0    0   0];
% C00*[x1^2; x2^2; x3^2]
C00 = scale_3* [0    0   0;
                0    0   0;
                0    0   0];
%> f(t,x) = A(mu,eta)x  + mu*M;
f_ls                    = @(t,y) A*y + M*mu + mu^2*A00*y + B00*(y.*y(index)) + C00*y.^2;
%
prob.odes_Fcn  = f_ls;

%> the period is around 6 second, so run the simulation for 100 cycles
prob.tspan  = [0 8000];
%>
prob.fs     = 20;

%> 
 prob.C     = C;
%> define  the function to get the accelaration in the observed coordinate
prob.ax_Fcn             = @(t,y) C*get_DF_DX(f_ls, y, [1 2 3])*f_ls(0,y);
%
prob.Wx_Fcn             = @(y)  -B;
prob.bx_Fcn             = @(y) C*get_DF_DX(f_ls, y, [1 2 3])*B;
prob.rx_Fcn             = @(x) par(4)+par(8)-1;
prob.Impact_map         = @(x) x - B*C*f_ls(0,x); 
%
prob.Post_impact_Fcn = @Shilnikov_homoclinic_chaos_post_impact;

%% define options
refine = 1;
prob.options = odeset('Events',@(t,y) ...
    Shilnikov_homoclinic_chaos_DetectingCross_events(t, y , prob),...
    'RelTol',1e-12,'AbsTol',1e-12,'Refine',refine);

prob.opsl_stuck=odeset(prob.options,'Events',@(t,y) ...
     Shilnikov_homoclinic_chaos_unstickEvent(t,y, prob));

end