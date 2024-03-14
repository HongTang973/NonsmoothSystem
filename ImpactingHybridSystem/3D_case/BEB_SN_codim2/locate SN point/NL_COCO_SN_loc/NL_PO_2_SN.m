%% > Shilnikov example to detect the PD point and show the codim-2 curve
clc
close all
clear 
original=pwd;
isUoB = get_platform();
%>
if isUoB; onedrive_path = 'C:\Users\ib20968\'; else; onedrive_path = 'F:\onedrive\'; end
cd([onedrive_path,'OneDrive - University of Bristol\Codes stall\coco_2020Mar22\coco'])
startup;
cd(original);
%
% Construct initial solution guess from already computed brute-force
% simulation: the initial condition
index 		=  [7;8]; %> the location of the flight velocity and the r
ds     		=  [1;0]; %> varying the parameter which controls the BEB

%% define the system
global const
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
const.A00   = A00;
const.B00   = B00;
const.C00   = C00;
%> this is chosen from case2 NL case, where there is a P1 cycle with chosen
%> parameter
% par = {'rho', 'omega', 'lambda','r0', 'c0', 'mu','eta'}
%% point 1
par_0 = [-0.1,0.2,-0.5,1.78192697901049,1.6,4.79596124049338,0.02,0.0955685096979647]';
%> IC
x0    = [0.131788909327915;0.117382181376020;0.00182991263983379];
x00   = [0; 0.0770926271121134; 0.0607559985948500];

%% point 2



%> define the autonoumous ode function with the input as: y,p,const,STRING
% Shilnikov_IHS_ode
%
% par    = [];
% modes  = {'free'   'free'};
% events = {'zero_v' 'impact'};
% resets = {'zero_v' 'bounce'};
par    = [];
modes  = {'free'};
events = {'impact'};
resets = {'bounce'};
%
f       = @(t, x) full_BEB_SN_IHS_3D_ode(x, par_0,'free');
% from the H(x) = 0 to the zero_v 
% [t1, x1] = ode45(f, [0 2.926150229406630], x00);
[t1, x1] = ode45(f, [0 6], x00);
% [t1, x1] = ode45(f, [0 2.9], x00);
% from the zero_v to the  H(x) = 0 
% [t2, x2] = ode45(f, [0 2.706], x0);
[t2, x2] = ode45(f, [0 2.70], x0);
figure
subplot(121)
plot(x1(:,1),x1(:,2))
subplot(122)
plot(x2(:,1),x2(:,2))
%
t0     = {t1 };
x0     = {x1 };
% t0     = {t1 t2};
% x0     = {x1 x2};

% Use anonymous, in-line definitions of the event and reset functions and
% their Jacobians.

% Shilnikov_IHS_events     = @(x,p,~) x(1);
% Shilnikov_IHS_events_dx  = @(x,p,~) [1 0 0];
% Shilnikov_IHS_events_dp  = @(x,p,~) zeros(1,8);
% Shilnikov_IHS_resets     = @(x,p,~) [x(1);x(2);0];
% Shilnikov_IHS_resets_dx  = @(x,p,~) [1 0 0; 0 1 0; 0 0 0];
% Shilnikov_IHS_resets_dp  = @(x,p,~) zeros(3,5);


%% initialize the tool and customized settings
prob=coco_prob();
% prob = coco_set(prob, 'ode', 'vectorized', true); % vectorized or not the Jacobian
prob = coco_set(prob,'cont','h0',0.01); % 0.0001 or 0.01 different results
prob = coco_set(prob,'cont','h_min',10^-2); %0.0001
prob = coco_set(prob,'cont','h_max',0.01); %1
prob = coco_set(prob,'cont','h_fac_min',0.1); % Minimum step size adaptation factor ,0.5 100
prob = coco_set(prob,'cont','h_fac_max',2); % Maximum step size adaptation factor, 2 35 (can do strange stuff if too low..like not escaping)
prob = coco_set(prob,'cont','almax', 65);   % Critical angle between successive tangent vectors (try changing this one), 35
prob = coco_set(prob, 'cont', 'NullItMX', 1);
prob = coco_set(prob, 'coll', 'NTST', 300); %mesh 30, 100, 10000
prob = coco_set(prob, 'coll', 'NCOL', 5); %mesh 5,7
prob = coco_set(prob, 'coll', 'TOL', 10^-2); %mesh
% prob = coco_set(prob, 'cont', 'PtMX', 50); % maximum number of continuation steps
prob = coco_set(prob, 'cont', 'ItMX', 5000); % maximum number of iteration - this is the one that impose the number of iterations
prob = coco_set(prob, 'cont', 'NPR', 1); %iteration
prob = coco_set(prob, 'cont', 'NAdapt', 2); % aadaptive changes are made to the orbit discretization after each successful step of continuation
prob = coco_set(prob,'corr','ItMX',50); %maximum number of iteration - for convergence (10 default)
prob = coco_set(prob,'corr','SubItMX',4); % damping in Newton criterion
prob = coco_set(prob,'corr','TOL',1.00E-006); %Tolerance of Newton correction
prob = coco_set(prob,'corr','ResTOL',1.00E-006); %Converge criterion of norm of the residium


% prob = ode_isol2hspo(prob, '', ...
%   {@Shilnikov_IHS_ode,  @Shilnikov_IHS_events,    @Shilnikov_IHS_resets}, ...
%   {@Shilnikov_IHS_DFDX, @Shilnikov_IHS_events_dx, @Shilnikov_IHS_resets_dx}, ...
%   {@Shilnikov_IHS_DFDP, @Shilnikov_IHS_events_dp, @Shilnikov_IHS_resets_dp}, ...
%   modes, events, resets, t0, x0, {'rho' 'ome' 'la','r0','c0', 'T','mu' 'eta'}, par_0);

prob = ode_isol2hspo(prob, '', ...
  {@full_BEB_SN_IHS_3D_ode,  @full_BEB_SN_IHS_3D_events,    @full_BEB_SN_IHS_3D_resets}, ...
  modes, events, resets, t0, x0, {'la_1' 'la_2' 'la_3','b2_0','b3', 'T','mu' 'eta'}, par_0);

% The @duff_add_IP slot function responds to the 'bddat' signal and stores
% the components of the initial end point on the first trajectory segment
% to the cell array returned by the coco entry-point function and stored to
% disk.

prob = coco_add_slot(prob, 'duff_bddat', @full_BEB_SN_IHS_3D_add_IP, [], 'bddat');
% prob = coco_set(prob, 'cont', 'h', 2, 'h_max', 10, 'PtMX', 200, 'NAdapt', 5);

fprintf('\n Run=''%s'': Continue primary family of two-segment periodic orbits.\n', ...
  'run1');

bd0 = coco(prob, 'run1', [], 1, 'mu', [-0.02 0.02]);

%% Continue along secondary branches of two-segment periodic orbits

% The encoding below is identical to the one above.

bd1  = coco_bd_read('run1');
labs = coco_bd_labs(bd1, 'SN');
for lab=labs
  prob = coco_prob();
  prob = ode_BP2hspo(prob, '', 'run1', lab);
  prob = coco_add_slot(prob, 'duff_bddat', @Shilnikov_IHS_add_IP, [], 'bddat');
  prob = coco_set(prob, 'cont', 'h', 2, 'h_max', 10, 'PtMX', 50, 'NAdapt', 5);
  
  fprintf(...
  '\n Run=''%s'': Continue secondary family of two-segment periodic orbits from point %d in run ''%s''.\n', ...
  sprintf('run2_%d',lab), lab, 'run1');

  coco(prob, sprintf('run2_%d',lab), [], 1, 'om', [0.5 1.5]);
end

%% Continuation of four-segment periodic orbits

% The encoding below uses solution data stored with a period-doubling
% bifurcation located along a secondary branch to restart continuation of a
% four-segment periodic orbit. The dimensional deficit is again 0,
% requiring the release of a single continuation parameter in order to
% obtain a one-dimensional solution manifold.

run = sprintf('run2_%d',labs(1));
bd2 = coco_bd_read(run);
labs = coco_bd_labs(bd2, 'PD');
for lab=labs([1 end])
  prob = coco_prob();
  prob = ode_PD2hspo(prob, '', run, lab);
  prob = coco_add_slot(prob, 'duff_bddat', @duff_add_IP, [], 'bddat');
  prob = coco_set(prob, 'cont', 'h', 2, 'h_max', 10, 'PtMX', 50, 'NAdapt', 5);
  
  fprintf(...
    '\n Run=''%s'': Continue family of period-doubled four-segment periodic orbits from point %d in run ''%s''.\n', ...
    sprintf('run3_%d', lab), lab, run);

  coco(prob, sprintf('run3_%d', lab), [], 1, 'om', [0.5 1.5]);
end

%% Continuation of period-doubling bifurcations

% The encoding below imposes additional constraints on the two-segment
% boundary-value problem corresponding to a period-doubling bifurcation,
% resulting in a dimensional deficit of -1. A one-dimensional solution
% manifold results by releasing 'om' and 'A' and allow these to vary during
% continuation.

bd1  = coco_bd_read('run1');
labs = coco_bd_labs(bd1, 'BP');
run  = sprintf('run2_%d',labs(1));
bd2  = coco_bd_read(run);
labs = coco_bd_labs(bd2, 'PD');
prob = coco_prob();
prob = coco_set(prob, 'cont', 'h', 2, 'h_max', 20, 'PtMX', 100, 'NAdapt', 1);
prob = coco_add_slot(prob, 'duff_bddat', @duff_add_IP, [], 'bddat');
prob = ode_PD2PD(prob, '', run, labs(1));

fprintf(...
  '\n Run=''%s'': Continue family of period-doubling bifurcations from point %d in run ''%s''.\n', ...
  'run4', labs(1), run);

coco(prob, 'run4', [], 1, {'om' 'A'}, [0.5 1.5]);

%% Continuation of saddle-node bifurcations

% The encoding below imposes additional constraints on the two-segment
% boundary-value problem corresponding to a saddle-node bifurcation,
% resulting in a dimensional deficit of -1. A one-dimensional solution
% manifold results by releasing 'om' and 'A' and allow these to vary during
% continuation.

bd1  = coco_bd_read('run1');
labs = coco_bd_labs(bd1, 'SN');
prob = coco_prob();
prob = coco_set(prob, 'cont', 'h', 2, 'h_max', 20, 'PtMX', 200, 'NAdapt', 5);
prob = coco_add_slot(prob, 'duff_bddat', @duff_add_IP, [], 'bddat');
prob = ode_SN2SN(prob, '', 'run1', labs(1));

fprintf(...
  '\n Run=''%s'': Continue family of saddle-node bifurcations from point %d in run ''%s''.\n', ...
  'run5', labs(1), 'run1');

coco(prob, 'run5', [], 1, {'om' 'A'}, [0.2 1.5]);

%% Graphical representation of results

% Panel (a)
bd1 = coco_bd_read('run1');

figure(1); clf; hold on; grid on; box on; view(160,25)

thm = struct('special', {{'EP', 'FP', 'BP', 'PD'}});
coco_plot_bd(thm, 'run1', 'om', 'A', 'X0')

labs = coco_bd_labs(bd1, 'BP');

for lab=labs
  coco_plot_bd(thm, sprintf( 'run2_%d',lab), 'om', 'A', 'X0')
end

run = sprintf('run2_%d',labs(1));
bd2 = coco_bd_read(run);
labs = coco_bd_labs(bd2, 'PD');

for lab=labs([1 end])
  coco_plot_bd(thm, sprintf('run3_%d',lab), 'om', 'A', 'X0')
end

thm = struct('special', {{'EP', 'FP'}}, ...
  'xlab', '\omega', 'zlab', 'x_1(0)');
coco_plot_bd(thm, 'run4', 'om', 'A', 'X0')
coco_plot_bd(thm, 'run5', 'om', 'A', 'X0')

hold off

% Plot data: panels (b)-(d)
labs  = [9 12 16];

seg1.sol.RO = {'LineStyle', '-', 'LineWidth', 2, 'Color', 'black', ...
  'Marker', '.', 'MarkerSize', 12};
seg2.sol.RO = {'LineStyle', '-', 'LineWidth', 2, 'Color', ...
  [0.6 0.6 0.6], 'Marker', '.', 'MarkerSize', 12};

for i=1:3
  figure(i+1); clf; hold on; box on; grid on
  coco_plot_sol(seg1, 'run1', labs(i), 'hspo.orb.bvp.seg1', 'x', 'x')
  coco_plot_sol(seg2, 'run1', labs(i), 'hspo.orb.bvp.seg2', 'x', 'x')
  axis([-inf inf -inf inf]); hold off
end


