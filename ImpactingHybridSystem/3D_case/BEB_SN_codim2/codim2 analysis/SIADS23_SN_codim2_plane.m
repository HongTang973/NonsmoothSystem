%> use the CO method to track the solution branch
close all
%> the initial point on the SN curve is from the file locating_SN_point.m
par_SN_init = [-0.1 + 0.2i	-0.1 - 0.2i	-0.5 1.78192697901049 	1.6	5.84686867166138 ];

%> the point hitted by the line continuation
% par_SN_init = [-0.1	0.2	-0.525014022367497	1.87753996111232	1.6	4.795364374429504];

prob.par    = par_SN_init;
indexes     = [3, 4, 6];
[t_,Jacob]= get_jacobian(@SIADS23_SN_zero_Fcns, par_SN_init, indexes) %> & done
% prob.keys   =    keys;
% prob.par2prob                   = @par2NForm_DummyVar;
prob.par2prob                   = @par2NForm_Lienard;
det_ = SIADS23_SN_zero_Fcns(par_SN_init)
event_val = SN_PD_MonitorFcns(prob,par_SN_init)
%  ------- new event detection functionality -------
prob.bifur_detct                = 0;
prob.Events.NumOfEvents         = 6;
prob.Events.MonitorFunctions    = @SN_PD_MonitorFcns;
prob.Events.tag                 = {'SN';'SN';'SN';'PD';'PD';'PD'};
% ---------------------------------------------------
prob.N_points   = 100;
prob.uplim_step = 4000;
prob.ind_p  =    4;
prob.ind_x  =    [3,6];
prob.dir    =    1;
prob.min_h  =    1e-4;
prob.co_dh  =    1e-1;
prob.p_span = [0.8*par_SN_init(4) 1.1*par_SN_init(4)];
prob.x_span = [-0.6 -0.4; 3 10];
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
%> the approximation for the (mu,eta) curve
% f_mu_eta =@(mu,eta) a0*mu + b0*eta - (a1*mu + b1*eta)^2/(4*c);
x_t_range_1 = [-0.02:0.0001:0.02];
x_t_range_2 = [-0.02:0.0001:0];
FIG1 = figure;
plot(u1(3,:),u1(4,:),'b-','linewidth', 1.5);
hold on
plot(u2(3,:),u2(4,:),'b-','linewidth', 1.5);
plot(par_SN_init(3),par_SN_init(4),'ro','MarkerSize',5,'MarkerFaceColor',[1 1 1])
xlabel('$\lambda_3$', 'interpreter','latex')
ylabel('$b_2$', 'interpreter','latex')

% legend([h1 h2 h3], 'Interpreter','latex')
% fimplicit(f_mu_eta,[-0.2 0.2 -0.2 0.2])
xlim([-0.55 -0.45])
% ylim([-0.2 0.1])
% grid on
adj_plot_theme_I(FIG1)
% exportgraphics(FIG1,'./Codim2_3D_case.pdf','ContentType','vector')
virtual_SN_mu  = u2(3,:)-par_SN_init(3);
virtual_SN_eta = u2(4,:)-par_SN_init(4);
admiss_SN_mu   = u1(3,:)-par_SN_init(3);
admiss_SN_eta  = u1(4,:)-par_SN_init(4);
FIG2 = figure;
h1 = plot(-virtual_SN_mu, virtual_SN_eta,'r-','linewidth', 1.5, 'displayname','Addmissible LCO');
hold on
h2 = plot(-admiss_SN_mu, admiss_SN_eta,'k--','linewidth', 1.5,'displayname','Virtual LCO');
h3 = plot(-x_t_range_1,t_line(x_t_range_1),'b-','linewidth', 0.8, 'displayname', '$\eta = t_2/t_1 \mu$');
plot(0,0,'ro','MarkerSize',5,'MarkerFaceColor',[1 1 1])
mu_crit     = 0.025014022367497;
eta_crit    = 0.095612982101833;
plot(mu_crit,eta_crit,'bd','MarkerSize',8,'MarkerFaceColor',[1 1 1])
xlabel('$\mu$', 'interpreter','latex')
ylabel('$\eta$', 'interpreter','latex')
legend([h1 h2 h3], 'Interpreter','latex','location','best')
grid on
adj_plot_theme_I(FIG2)
% exportgraphics(FIG2,'./Codim2_3D_case_mu_eta.pdf','ContentType','vector')


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
