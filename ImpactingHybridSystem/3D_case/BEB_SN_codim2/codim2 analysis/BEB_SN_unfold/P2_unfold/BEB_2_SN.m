%% use the continuation to track the stable and the unstable LCO from BEB
% to SN point -> then save the data as BEB_2_SN_LCO_continuation.mat for
% the BEB_SN_P2_unfold.m to use


clc
clear
close all

%>  start from the middle point
P2 = [-0.1,0.2,-0.5,1.78192697901049,1.6,4.79596124049338,0.025,0.0955685096979647];
%> slight variation from the SN 
% P2(7) = 0.5*P2(7);
mu_crit      = P2(7);
eta_crit     = P2(8);
T_MP         = 3.379162984466553;
par_MP_Nform = [-0.1, 0.2, -0.5, 1.781926979010490, 1.6, T_MP, 0.5*mu_crit, eta_crit];

%% initial guess for the two roots
% t_bt      = 0;
% t_up      = 10;
% t_range   = t_bt:1/1e4:t_up;
% 
% par_unfold = par_MP_Nform;
% [A,B,C,R,T_2_det] = par2NForm_DummyVar(par_unfold);
% det_range = T_2_det(t_range );
% figure
% plot(t_range, det_range)
% hold on 
% plot(t_range, 0*t_range, 'r--')
% plot(haxes1,par_unfold(7),par_unfold(8),'gs')
US_T        = 3.38;
S_T         = 6.674;
par_MP_US = par_MP_Nform; par_MP_US(6)=US_T;
par_MP_S  = par_MP_Nform; par_MP_S(6) =S_T;


% prob.keys   =    keys;
prob.par2prob   = @par2NForm_DummyVar;
%  ------- new event detection functionality -------
prob.bifur_detct                = 1;
prob.Events.NumOfEvents         = 6;
prob.Events.MonitorFunctions    = @SN_PD_MonitorFcns;
% prob.Events.MonitorFunctions    = @SIAMDS23_SN_codim_MonitorFcns;
prob.Events.tag                 = {'SN';'SN';'SN';'PD';'PD';'PD'};
% ---------------------------------------------------
prob.N_points       = 100;
prob.uplim_step     = 6000;
prob.ind_p          =    7;
prob.ind_x          =    6;
prob.dir            =    1;
prob.co_dh          =    1e-2;
prob.min_h          =    1e-3;
prob.p_span         = [0 0.03];
prob.x_span         = [0 10];
zero_det            = @(par) [1,0]*SIADS23_SN_zero_Fcns_Dummy(par);
prob.ZeroFunctions  = zero_det;
%

% F = SIADS23_SN_zero_Fcns_Dummy( prob.par)
% % F = SIADS23_SN_zero_Fcns(par_SP_Nform)
% indexes     = [6, 7, 8];
% [t_,Jacob]  = get_jacobian(@SIADS23_SN_zero_Fcns_Dummy, prob.par, indexes) %> & done
% t_(3)/t_(2)
%% track the stable branch
% par_MP_S = [-0.1;0.2;-0.5;1.78192697901049;1.6;2.91571041661273;0.1;0.0955685096979647];
prob.par        = par_MP_S; 
tic
prob_1 = prob;
prob_2 = prob;
prob_2.dir    =    -1;
[u1,iter1,lab1,Monitor_fval_1] = codim1_PC(prob_1);
[u2,iter2,lab2,Monitor_fval_2] = codim1_PC(prob_2);
toc
SN_locus_1    = strcmp(lab1,'SN');
SN_locus_2    = strcmp(lab2,'SN');

PD_locus_1    = strcmp(lab1,'PD');
PD_locus_2    = strcmp(lab2,'PD');

S_locus_1    = strcmp(lab1,'S');
S_locus_2    = strcmp(lab2,'S');

US_locus_1    = strcmp(lab1,'US');
US_locus_2    = strcmp(lab2,'US');
%>
US_u          = [u1(:,US_locus_1), u2(:,US_locus_2)];
S_u           = [u1(:,S_locus_1), u2(:,S_locus_2)];
%> resort the sequence
[~,index_US]       = sort(US_u(7,:),'ascend');
US_u               = US_u(:,index_US);
[~,index_S]       = sort(S_u(7,:),'ascend');
S_u               = S_u(:,index_S);
%% track the unstable branch
% prob.par        = par_MP_US; 
% tic
% prob_1 = prob;
% prob_2 = prob;
% prob_2.dir    =    -1;
% [US_u1,US_iter1,US_lab1,US_Monitor_fval_1] = codim1_PC(prob_1);
% [US_u2,US_iter2,US_lab2,US_Monitor_fval_2] = codim1_PC(prob_2);
% toc

FIG1 = figure;
haxes1 = axes;
h1 = plot(S_u(7,:),S_u(6,:),'r-','linewidth', 1.5, 'displayname','Stable LCO');
hold on
h2 = plot(US_u(7,:),US_u(6,:),'k--','linewidth', 1.5,'displayname','Unstable LCO');

