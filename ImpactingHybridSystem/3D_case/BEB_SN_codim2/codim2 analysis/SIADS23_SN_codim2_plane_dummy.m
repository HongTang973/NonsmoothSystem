%% > use the CO method to track the solution branch
clc
close all
par_SN_init     = [-0.1, 0.2, -0.5, 1.781926979010490, 1.6, 5.846868671661380];
prob.par        = [par_SN_init,0,0]; 
prob.par2prob   = @par2NForm_DummyVar;
F = SIADS23_SN_zero_Fcns_Dummy( prob.par)
% F = SIADS23_SN_zero_Fcns(par_SP_Nform)
indexes     = [6, 7, 8];
[t_,Jacob]  = get_jacobian(@SIADS23_SN_zero_Fcns_Dummy, prob.par, indexes) %> & done
t_(3)/t_(2)
% prob.keys   =    keys;
%%  ------- new event detection functionality -------
prob.bifur_detct                = 1;
prob.Events.NumOfEvents         = 6;
% prob.Events.MonitorFunctions    = @SN_PD_MonitorFcns;
prob.Events.MonitorFunctions    = @SIAMDS23_SN_codim_MonitorFcns;
prob.Events.tag                 = {'PD';'PD';'PD';'CR';'CR';'CR'};
% ---------------------------------------------------
prob.N_points   = 100;
prob.uplim_step = 2000;
prob.ind_p  =    7;
prob.ind_x  =    [6,8];
prob.dir    =    1;
prob.co_dh  =    1e-2;
prob.min_h  =    1e-4;
prob.p_span = [-0.06 0.06];
prob.x_span = [4 10; -0.6 0.4];
prob.ZeroFunctions = @SIADS23_SN_zero_Fcns_Dummy;

tic
prob_1 = prob;
prob_2 = prob;
prob_2.dir    =    -1;
[u1,iter1,lab1,Monitor_fval_1] = codim1_PC(prob_1);
[u2,iter2,lab2,Monitor_fval_2] = codim1_PC(prob_2);
toc
% 
CR_locus_1    = strcmp(lab1,'CR');
CR_locus_2    = strcmp(lab2,'CR');

CR_points           = [u1(:,CR_locus_1)'; u2(:,CR_locus_2)'];
%> 
t_line = @(x) t_(3)/t_(2)*x;
x_t_range_1 = [-0.04:0.001:0.04];
x_t_range_2 = [-0.02:0.0001:0];
FIG1 = figure;
haxes1 = axes;
h1 = plot(u1(7,:),u1(8,:),'r-','linewidth', 1.5, 'displayname','Addmissible LCO');
hold on
h2 = plot(u2(7,:),u2(8,:),'k--','linewidth', 1.5,'displayname','Virtual LCO');
h3 = plot(x_t_range_1,t_line(x_t_range_1),'b-','linewidth', 0.8, 'displayname', '$\eta = t_2/t_1 \mu$');
plot(0,0,'ro','MarkerSize',5,'MarkerFaceColor',[1 1 1])
% plot(x_t_range_2,t_line(x_t_range_2),'b--','linewidth', 1.2)
xlabel('$\mu$', 'interpreter','latex')
ylabel('$\eta$', 'interpreter','latex')
legend([h1 h2 h3], 'Interpreter','latex','location','best')
xlim([-0.06 0.06])
ylim([-0.15 0.15])
grid on
adj_plot_theme_I(FIG1)
% exportgraphics(FIG1,'./Codim2_3D_case.pdf','ContentType','vector')
% figure
% plot(u1(7,:), u1(6,:))
% hold on
% plot(u2(7,:), u2(6,:))

virtual_SN_mu  = u2(7,:);
virtual_SN_eta = u2(8,:);
admiss_SN_mu   = u1(7,:);
admiss_SN_eta  = u1(8,:);

FIG2 = figure;
plot(virtual_SN_mu, virtual_SN_eta)
hold on
plot(admiss_SN_mu, admiss_SN_eta)



%>
t_bt      = 0;
t_up      = 10;
t_range   = t_bt:1/1e4:t_up;


% [A,B,C,R, T_2_det] = par2NForm(par_SN_init);
% det_range = T_2_det(t_range );
% figure
% plot(t_range, det_range)
% hold on 
% plot(t_range, 0*t_range, 'r--')

par_sw = [par_SN_init, 0.01, -0.1];
[A,B,C,R,T_2_det] = par2NForm_DummyVar(par_sw);
det_range = T_2_det(t_range );
figure
plot(t_range, det_range)
hold on 
plot(t_range, 0*t_range, 'r--')
plot(haxes1,par_sw(7),par_sw(8),'ro')

par_ne = [par_SN_init, 0.01, 0];
[A,B,C,R,T_2_det] = par2NForm_DummyVar(par_ne);
det_range = T_2_det(t_range );
figure
plot(t_range, det_range)
hold on 
plot(t_range, 0*t_range, 'r--')
plot(haxes1,par_ne(7),par_ne(8),'ro')
%>
t_bt      = 0;
t_up      = 10;
t_range   = t_bt:1/1e4:t_up;


% [A,B,C,R, T_2_det] = par2NForm(par_SN_init);
% det_range = T_2_det(t_range );
% figure
% plot(t_range, det_range)
% hold on 
% plot(t_range, 0*t_range, 'r--')

par_sn = [par_SN_init, 0.02, 0.1];
[A,B,C,R,T_2_det] = par2NForm_DummyVar(par_sn);
det_range = T_2_det(t_range );
figure
plot(t_range, det_range)
hold on 
plot(t_range, 0*t_range, 'r--')
plot(haxes1,par_sn(7),par_sn(8),'gs')

% mu_crit     = 0.025014022367497;
% eta_crit    = 0.095612982101833;
% par_unfold = [par_SN_init, 0.02, 0.1];
par_unfold = CR_points(2,:);
[A,B,C,R,T_2_det] = par2NForm_DummyVar(par_unfold);
det_range = T_2_det(t_range );
figure
plot(t_range, det_range)
hold on 
plot(t_range, 0*t_range, 'r--')
plot(haxes1,par_unfold(7),par_unfold(8),'gs')
F = SIADS23_SN_zero_Fcns_Dummy( par_unfold)


