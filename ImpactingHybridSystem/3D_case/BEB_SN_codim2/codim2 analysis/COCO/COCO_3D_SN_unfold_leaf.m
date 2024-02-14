%% build-up of the coco in the non-smooth impacting hybrid system
% clear
% clc
close all
%% -------------  get the coco's directory ----
% original=pwd;
% valve_project_startup;
% cd(original);

%% initialize the tool and customized settings
prob=coco_prob();
% prob = coco_set(prob, 'ode', 'vectorized', true); % vectorized or not the Jacobian
prob = coco_set(prob,'cont','h0',0.001); % 0.0001 or 0.01 different results
prob = coco_set(prob,'cont','h_min',10^-4); %0.0001
prob = coco_set(prob,'cont','h_max',0.02); %1
prob = coco_set(prob,'cont','h_fac_min',0.1); % Minimum step size adaptation factor ,0.5 100
prob = coco_set(prob,'cont','h_fac_max',2); % Maximum step size adaptation factor, 2 35 (can do strange stuff if too low..like not escaping)
prob = coco_set(prob,'cont','almax', 45);   % Critical angle between successive tangent vectors (try changing this one), 35
prob = coco_set(prob, 'cont', 'NullItMX', 1);
prob = coco_set(prob, 'coll', 'NTST', 30); %mesh 30, 100, 10000
prob = coco_set(prob, 'coll', 'NCOL', 5); %mesh 5,7
prob = coco_set(prob, 'coll', 'TOL', 10^-2); %mesh
% prob = coco_set(prob, 'cont', 'PtMX', 50); % maximum number of continuation steps
prob = coco_set(prob, 'cont', 'ItMX', 500); % maximum number of iteration - this is the one that impose the number of iterations
prob = coco_set(prob, 'cont', 'NPR', 1); %iteration
prob = coco_set(prob, 'cont', 'NAdapt', 1); % aadaptive changes are made to the orbit discretization after each successful step of continuation
prob = coco_set(prob,'corr','ItMX',50); %maximum number of iteration - for convergence (10 default)
prob = coco_set(prob,'corr','SubItMX',4); % damping in Newton criterion
prob = coco_set(prob,'corr','TOL',1.00E-009); %Tolerance of Newton correction
prob = coco_set(prob,'corr','ResTOL',1.00E-009); %Converge criterion of norm of the residium

% ------------- step 1: the zero problem ---
% Determine where your m-file's folder is.
% folder = fileparts(which(Obeservation_Funs_Airfoil.m));
% folder = 'F:\onedrive\OneDrive - University of Bristol\Codes stall\3DOF_AIRFOIL_MODEL\Bifurcation_Analysis\BEB\PO\';
% % Add that folder plus all subfolders to the path.
% addpath(genpath(folder));

%% continuation of the equlibrium branch
%  ----------   adding zero functions    ----------------------- %
% the command is as :
% ### coco_add_func(proc, fid, fhan, data, 'zero', 'uidx', uidx, 'u0', u0)
% ---fid: contains a function identifier
% ---fhan: contain a function handle
% ---data: contain the function data structure
% ---uidx: contain the index set
% ---u0: contain the vector u0
% ---'zero': function type

% ### coco_read_solution(fid, runid, labelID)
%> the parameter is chosen on the co-dim 2 curve: via the SN_codim2_plane
mu_crit     = 0.025014022367497;
eta_crit    = 0.095612982101833;
T_crit      = 4.795364374429504;
par_SP_Nform        = [-0.1, 0.2, -0.5, 1.781926979010490, 1.6, T_crit, mu_crit, eta_crit];
T_MP         = 3.379162984466553;
%> the initial point is selected away from the fold point
par_MP_Nform = [-0.1, 0.2, -0.5, 1.781926979010490, 1.6, T_MP, 0.5*mu_crit, eta_crit];
F = SIADS23_SN_zero_Fcns_Dummy(par_MP_Nform)
F = SIADS23_SN_zero_Fcns_Dummy(par_SP_Nform)

%> the Floquet multiplier
[A,B,C,R,T_2_det]   = par2NForm_DummyVar(par_SP_Nform);
x_00                = IC_generator(par_SP_Nform(6),R,A,C,1);
L0                  = C*A*inv(R)*x_00;
[Mono_p,Salt_p]     = IC2Floque_Multipliers(par_SP_Nform(6),x_00,R,A,C);
T_2_det(par_SP_Nform(6))
%> solve the pseudo equilibrium
hat_y = pinv([C;C*A;(eye(length(C)) - B*C*A/(C*A*B))])*[1;0;0*C']
%> From the parameter point near the SN 
par                 =  par_MP_Nform;
keys                = {'lambda_1', 'lambda_2', 'lambda_3','b2','b3','T', 'mu','eta'};
prob.par2prob       = @par2NForm_DummyVar;
%> pass the initial condition to the COCO
u0                  = par;
%> initialize the prob with customized settings



%% > embed the valve problem into the coco continuation problem

%> construct the problem from COCO built-in function
prob = coco_add_func(prob, 'Det_func', @COCO_3D_LCO_IHS, [], 'zero', 'u0', u0);

p_span   = [0.0 0.4];
X_span   = [0 10];
% define continuation parameters
PNM             = keys;
PIdx1           = [1;2;3;4;5;6;7;8];
prob            = coco_add_pars(prob, '' ,PIdx1, PNM(PIdx1),'inactive');
% bd0 = coco(prob, '3D_case_mu_T', [], 1, {'mu','T'},  {p_span,X_span});

%>
% figure; clf
% haxes1 = axes;
% thm = struct('special', {{'BP','FP'}});
% coco_plot_bd(thm, '3D_case_mu_T', 'mu', 'T')
% hold on

%> observe the evolution of the fixed point
mu_list_bd0             = coco_bd_col(bd0, {'mu'});
T_list_bd0              = coco_bd_col(bd0, {'T'});
FPlab_bd0               = coco_bd_labs(bd0, 'FP');
lab_list_bd0            = coco_bd_col(bd0, 'LAB');
bull_ismember_bd0       = ismember(lab_list_bd0,FPlab_bd0);
FP_index_bd0            = find(bull_ismember_bd0);
%>
n_num           = length(mu_list_bd0);
par_matrix      = repmat(par,n_num,1);
par_matrix(:,6) = T_list_bd0;
par_matrix(:,7) = mu_list_bd0;
%> get the fixed point 
[A,B,C,R] = par2NForm_DummyVar(par_SP_Nform);
% par_matrix(31,:); %> where the irregular point emerges
 x_00       = IC_generator(par_SP_Nform(6),R,A,C,1);
%  u_p        = x_00(2:3);
[Mono_p,Salt_p]     = IC2Floque_Multipliers(par_SP_Nform(6),x_00,R,A,C);
x_0    = R\x_00;
u_p    = x_0(2:3);
% GM_p   = max(abs(Salt_p))
% v = [-0.659548604166351;-0.751661917847524];
% w = [-0.561944587291412;-0.837304933094644];
hat_u  = [];
proj_u = [];
norm_u = [];
GM     = [];
for i =1 : n_num
    [A,B,C,R] = par2NForm_DummyVar(par_matrix(i,:));
    x_00                = IC_generator(par_matrix(i,6),R,A,C,1);
    [Mono_p,Salt_p]     = IC2Floque_Multipliers(par_matrix(i,6),x_00,R,A,C);
    GM                  = [GM, max(abs(Salt_p))];
%     hat_u       = [hat_u, x_00(2:3)];
%     norm_u      = [norm_u,norm(x_00(2:3) - u_p)];
%     proj_u      = [proj_u,w'*(x_00(2:3) - u_p)];
    
    x_0    = R\x_00;
    hat_u       = [hat_u, x_0(2:3)];
    norm_u      = [norm_u,norm(x_0(2:3) - u_p)];
    proj_u      = [proj_u,w'*(x_0(2:3) - u_p)];
end
%>  plot the distribution of the fixed points on the Sigma
figure;clf
plot(hat_u(1,:), hat_u(2,:),'b.')
hold on
plot(hat_u(1,GM>1+2e-5), hat_u(2,GM>1+2e-5),'r.')
plot(u_p(1), u_p(2), 'ro')
plot([u_p(1),u_p(1)+0.02*v(1)],[u_p(2),u_p(2)+0.02*v(2)],'k-' )
plot([u_p(1),u_p(1)-0.02*v(1)],[u_p(2),u_p(2)-0.02*v(2)],'k-' )

%> 
mu_at_FP_bd0            = real(mu_list_bd0(FP_index_bd0));
norm_at_FP_bd0          = real(norm_u(FP_index_bd0));

% figure;clf
% plot(mu_list_bd0(GM>1+2e-5), norm_u(GM>1+2e-5), 'r.')
% hold on
% plot(mu_list_bd0(GM<=1+2e-5), norm_u(GM<=1+2e-5), 'b.')
% plot(mu_at_FP_bd0, norm_at_FP_bd0, 'ro')
% 1.8463
% a= -0.200899639382884;
% c = 1.2208;
% z_curve_para = @(mu) -sqrt(-a/c*( mu_crit -mu));

%> added on 16th/Dec 
d      = 0;
% a1 =0;
c_     = @(mu) sqrt(3*a1*d*(mu_crit - mu) + c^2);
on_off = 1;
z_curve_para_1 = @(mu) sqrt(-a0./c_(mu).*(mu_crit - mu) + on_off*a1^2/4./c_(mu).^2.*(mu_crit - mu).^2) - on_off*a1/2./c_(mu).*(mu_crit - mu);
z_curve_para_2 = @(mu) -sqrt(-a0./c_(mu).*(mu_crit - mu) + on_off*a1^2/4./c_(mu).^2.*(mu_crit - mu).^2) - on_off*a1/2./c_(mu).*(mu_crit - mu);


FIG1 = figure;clf
plot(mu_list_bd0(GM>1+2e-5), proj_u(GM>1+2e-5), 'r.')
hold on
plot(mu_list_bd0(GM<=1+2e-5), proj_u(GM<=1+2e-5), 'b.')
plot(mu_at_FP_bd0, norm_at_FP_bd0, 'ro')
% plot(mu_list_bd0, z_curve_para(mu_list_bd0), 'k--')
% plot(mu_list_bd0, -z_curve_para(mu_list_bd0), 'k--')
plot(mu_list_bd0, z_curve_para_1(mu_list_bd0), 'g--')
plot(mu_list_bd0, z_curve_para_2(mu_list_bd0), 'k-')
xlabel('$\mu$','Interpreter','latex')
ylabel('$w^{\top}(\hat{u} - u_p)$','Interpreter','latex')
adj_plot_theme_I(FIG1)


%>
% figure
% plot(mu_list_bd0,   GM)
% hold on
% plot(mu_crit,       GM_p, 'ro')
% 
% figure
% plot(mu_list_bd0, proj_u)

%> 

FIG2 = figure;
L1   = v(2);
h1 = plot(mu_list_bd0(GM>1+2e-5), L0*mu_list_bd0(GM>1+2e-5) + L1*mu_list_bd0(GM>1+2e-5).*proj_u(GM>1+2e-5),'b--','LineWidth',1.4,'displayname','Unstable LCO');
hold on
h2 = plot(mu_list_bd0(GM<=1+2e-5), L0*mu_list_bd0(GM<=1+2e-5) + L1*mu_list_bd0(GM<=1+2e-5).*proj_u(GM<=1+2e-5),'b-','LineWidth',1.4,'displayname','Stable LCO');
h3 = plot(mu_list_bd0, L0*mu_list_bd0 + L1*mu_list_bd0.*z_curve_para_1(mu_list_bd0), 'r--','LineWidth',1.4,'displayname','Unstable LCO');
h4 = plot(mu_list_bd0, L0*mu_list_bd0 + L1*mu_list_bd0.*z_curve_para_2(mu_list_bd0), 'r-','LineWidth',1.4,'displayname','Stable LCO');
h5 = plot(mu_crit, L0*mu_crit + L1*mu_crit.*z_curve_para_2(mu_crit), 'kd','MarkerSize',5,'MarkerFaceColor',[1 1 1],'displayname','SN');
h6 = plot(mu_list_bd0, 0*mu_list_bd0, 'g-','LineWidth',2,'displayname','Stable PE');
h7 = plot(-mu_list_bd0, 0*mu_list_bd0, 'y-','LineWidth',2,'displayname','Stable AE');
h8 = plot(0,0,'ko','MarkerSize',5,'MarkerFaceColor',[1 1 1],'displayname','BEB');
legend([h1 h2 h3 h4 h5 h6 h7,h8],'location','best')
grid on
xlim([-0.005 0.03])
ylim([-2.5 0.1]*1e-3)
xlabel('$\mu$','Interpreter','latex')
ylabel('$\mathcal{A}$','Interpreter','latex')
adj_plot_theme_I(FIG2)
% exportgraphics(FIG2,'./Codim2_3D_unfold_OP.pdf','ContentType','vector')
FIG3 = figure;
h1 = plot(mu_list_bd0(GM>1+2e-5), hat_u(1,GM>1+2e-5)-u_p(1,:), 'b--','LineWidth',1.4,'displayname','Unstable LCO');
hold on
h2 = plot(mu_list_bd0(GM<=1+2e-5), hat_u(1,GM<=1+2e-5)-u_p(1,:), 'b-','LineWidth',1.4,'displayname','Stable LCO');

h3 = plot(mu_list_bd0,  v(1)*z_curve_para_1(mu_list_bd0), 'r--','LineWidth',1.4,'displayname','Unstable LCO');
h4 = plot(mu_list_bd0, v(1)*z_curve_para_2(mu_list_bd0), 'r-','LineWidth',1.4,'displayname','Stable LCO');

%>
h5 = plot(mu_crit, 0, 'ko','MarkerSize',5,'MarkerFaceColor',[1 1 1],'displayname','SN');

xlabel('$\mu$','Interpreter','latex')
ylabel('$\delta(\hat{\mathcal{A}})$','Interpreter','latex')
legend([h1 h2 h3 h4 h5],'location','best')
grid on
adj_plot_theme_I(FIG3)
% exportgraphics(FIG3,'./Codim2_3D_unfold_scaled.pdf','ContentType','vector')