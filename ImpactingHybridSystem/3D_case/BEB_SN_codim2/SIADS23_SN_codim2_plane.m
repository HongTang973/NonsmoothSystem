%> use the CO method to track the solution branch

%> the initial point on the SN curve is from the file locating_SN_point.m
par_SN_init = [-0.1 + 0.2i	-0.1 - 0.2i	-0.5 1.78192697901049 	1.6	5.84686867166138 ];
prob.par    = par_SN_init;
indexes     = [3, 4, 6];
[t_,Jacob]= get_jacobian(@SIADS23_SN_zero_Fcns, par_SN_init, indexes) %> & done
% prob.keys   =    keys;
prob.par2prob                   = @par2NForm_DummyVar;
% prob.par2prob                   = @par2NForm_Lienard;

%  ------- new event detection functionality -------
prob.bifur_detct                = 0;
prob.Events.NumOfEvents         = 6;
prob.Events.MonitorFunctions    = @SN_PD_MonitorFcns;
prob.Events.tag                 = {'SN';'SN';'SN';'PD';'PD';'PD'};
% ---------------------------------------------------
prob.uplim_step = 2000;
prob.ind_p  =    4;
prob.ind_x  =    [3,6];
prob.dir    =    1;
prob.co_dh  =    1e-1;
prob.p_span = [0.9*par_SN_init(4) 1.05*par_SN_init(4)];
prob.x_span = [-0.6 -0.4; 4 9];
prob.ZeroFunctions = @SIADS23_SN_zero_Fcns;

tic
prob_1 = prob;
prob_2 = prob;
prob_2.dir    =    -1;
% prob_1.p_span = [1e-2 1];
% prob_1.co_dh  = 1e-2;
[u1,iter1] = codim1_PC(prob_1);
[u2,iter2] = codim1_PC(prob_2);
toc
%> 
t_line = @(x) t_(2)/t_(1)*x;
x_t_range_1 = [-0.02:0.0001:0.02];
x_t_range_2 = [-0.02:0.0001:0];
FIG1 = figure;
h1 = plot(u1(3,:)-par_SN_init(3),u1(4,:)-par_SN_init(4),'r-','linewidth', 1.5, 'displayname','Addmissible LCO');
hold on
h2 = plot(u2(3,:)-par_SN_init(3),u2(4,:)-par_SN_init(4),'k--','linewidth', 1.5,'displayname','Virtual LCO');
h3 = plot(x_t_range_1,t_line(x_t_range_1),'b-','linewidth', 0.8, 'displayname', '$\eta = t_2/t_1 \mu$');
plot(0,0,'ro','MarkerSize',5,'MarkerFaceColor',[1 1 1])
% plot(x_t_range_2,t_line(x_t_range_2),'b--','linewidth', 1.2)
xlabel('$\mu$', 'interpreter','latex')
ylabel('$\eta$', 'interpreter','latex')
legend([h1 h2 h3], 'Interpreter','latex')
xlim([-0.05 0.05])
ylim([-0.2 0.1])
grid on
adj_plot_theme_I(FIG1)
% exportgraphics(FIG1,'./Codim2_3D_case.pdf','ContentType','vector')
virtual_SN_mu  = u2(3,:)-par_SN_init(3);
virtual_SN_eta = u2(4,:)-par_SN_init(4);
admiss_SN_mu   = u1(3,:)-par_SN_init(3);
admiss_SN_eta  = u1(4,:)-par_SN_init(4);
% figure
% plot(virtual_SN_mu, virtual_SN_eta)
% hold on
% plot(admiss_SN_mu, admiss_SN_eta)
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

par_ne = [par_SN_init, 0.01, 0];
[A,B,C,R,T_2_det] = par2NForm_DummyVar(par_ne);
det_range = T_2_det(t_range );
figure
plot(t_range, det_range)
hold on 
plot(t_range, 0*t_range, 'r--')


%
