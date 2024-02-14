%> before run: initialize the IHS analysis toolbox
%> run the IHS_startup.m

%> run the general toolbox to get the functionality 

%> check the root of the det function of given parameters in 3D case
clc
close all
%> the parameter set in the paper SIADS23
par                 = [-0.1000 + 0.2000i, -0.1000 - 0.2000i, -0.5 , 1.8,  1.6, 10];
%> temporary par to be used for unfolding the structure of (z,mu)
par                 =  [-0.1, 0.2, -0.492102081657598, 1.744171587227222, 1.6, 6.177398495078845];

par             = [-0.1, 0.2, -0.5, 1.781926979010490, 1.6, T_crit, 0.5*mu_crit, eta_crit];
%> codim-2 point case 0: 1.781916719621720
% par                   = [-0.1000 + 0.2000i, -0.1000 - 0.2000i, -0.5 , 1.781916719621720+0.00006,  1.6  , 5.847542373676547 ];

% par                   = [-0.1000 + 0.2000i, -0.1000 - 0.2000i, -0.5 , 1.781916719621720+0.000010259409,  1.6  , 5.847542373676547 ];
% par = u1(:,end);
%> codim-2 point case 1:
% par                 = [-0.1000 + 0.2000i, -0.1000 - 0.2000i, -0.5 , 1.8,  1.591103794399196  , 5.698285027182902 ];

%> codim-2 point case 2:
% par                 = [-0.1000 + 0.2000i, -0.1000 - 0.2000i, -0.504103038128386 , 1.8,  1.6  , 5.676046546984852 ];

keys                  = {'lambda_1', 'lambda_2', 'lambda_3','b2','b3','T'};
% [A,B,C,R,T_2_det] = par2NForm(par);
 [A,B,C,R,T_2_det] = par2NForm_full(par);
%> evoke the det condition for the existence of 3D case
% det_condition_for_3D_general_OB_canonical_form;
% fprintf('The analytic formula has been verified! \n')
%> 
t_bt      = 1;
t_up      = 8;
t_range   = t_bt:1/1e4:t_up;

det_range = T_2_det(t_range);

figure
plot(t_range, det_range)

[roots, fvals] =get_roots_Det(A,R,C, [0 t_up],1)
abs(roots(2) - roots(1))

% figure(1)
% hold on 
% plot(t_range, det_range, 'b--', 'linewidth', 1.2)
% ylim([-1 1]*10^(-9))
% xlim([t_bt t_up])



