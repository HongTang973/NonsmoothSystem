clc
close all
%>  fix the parameter eta and  track the solution branch 
% mu_crit     = admiss_SN_mu(267);
% eta_crit    = admiss_SN_eta(267);
% T_crit      = u1(6,267);
mu_crit     = 0.025014022367497;
eta_crit    = 0.095612982101833;
T_crit      = 4.795364374429504;

par_SP_Nform        = [-0.1, 0.2, -0.5, 1.781926979010490, 1.6, T_crit, mu_crit, eta_crit];
indexes             = [6, 7, 8];
[t_,Jacob]          = get_jacobian(@SIADS23_SN_zero_Fcns_Dummy, par_SP_Nform, indexes)
[A,B,C,R,T_2_det]   = par2NForm_DummyVar(par_SP_Nform);
T_2_det(par_SP_Nform(6))

% F_zero(par_SP_Nform)


%>
% det_range = T_2_det(t_range );
% figure
% plot(t_range, det_range)
% hold on 
% plot(t_range, 0*t_range, 'r--')
%> continuation from the mid point 
T_MP         = 3.379162984466553;
par_MP_Nform = [-0.1, 0.2, -0.5, 1.781926979010490, 1.6, T_MP, 0.5*mu_crit, eta_crit];
[A,B,C,R,T_2_det]   = par2NForm_DummyVar(par_MP_Nform);
T_2_det(par_MP_Nform(6))
SIADS23_SN_zero_Fcns_Dummy(par_MP_Nform)
%
prob.par        = par_MP_Nform;
prob.uplim_step = 5000;
prob.ind_p  =    8;
prob.ind_x  =    6;
prob.dir    =    1;
prob.co_dh  =    1e-1;
prob.p_span =    [-0.1 2];
prob.x_span =    [1 10];
zero_det    = @(par) [1,0]*SIADS23_SN_zero_Fcns_Dummy(par);
prob.ZeroFunctions  = zero_det;
prob.par2prob       = @par2NForm_DummyVar;
%>
prob.bifur_detct                = 1;
prob.Events.NumOfEvents         = 6;
prob.Events.MonitorFunctions    = @SN_PD_MonitorFcns;
prob.Events.tag                 = {'SN';'SN';'SN';'PD';'PD';'PD'};


%%>

tic
prob_1 = prob;
prob_2 = prob;
prob_2.dir    =    -1;
% prob_1.p_span = [1e-2 1];
% prob_1.co_dh  = 1e-2;
[u1,iter1,lab1] = codim1_PC(prob_1);
[u2,iter2,lab2] = codim1_PC(prob_2);
toc

SN_locus_1    = strcmp(lab1,'SN');
SN_locus_2    = strcmp(lab2,'SN');

PD_locus_1    = strcmp(lab1,'PD');
PD_locus_2    = strcmp(lab2,'PD');
%% >
figure;
plot(u1(8,:), u1(6,:))
hold on
plot(u2(8,:), u2(6,:))
plot(u1(8,SN_locus_1), u1(6,SN_locus_1),'kd','MarkerFaceColor',[1 1 1])
plot(u2(8,SN_locus_2), u2(6,SN_locus_2),'kd','MarkerFaceColor',[1 1 1])
plot(u1(8,PD_locus_1), u1(6,PD_locus_1),'bo','MarkerFaceColor',[1 1 1])
plot(u2(8,PD_locus_2), u2(6,PD_locus_2),'bo','MarkerFaceColor',[1 1 1])
xlabel('r')
ylabel('T')

% figure;
% plot(u1(5,:), u1(6,:))
% hold on
% plot(u2(5,:), u2(6,:))
% plot(u1(5,SN_locus_1), u1(6,SN_locus_1),'kd','MarkerFaceColor',[1 1 1])
% plot(u2(5,SN_locus_2), u2(6,SN_locus_2),'kd','MarkerFaceColor',[1 1 1])
% plot(u1(5,PD_locus_1), u1(6,PD_locus_1),'bo','MarkerFaceColor',[1 1 1])
% plot(u2(5,PD_locus_2), u2(6,PD_locus_2),'bo','MarkerFaceColor',[1 1 1])
% xlabel('c')
% ylabel('T')
%>
