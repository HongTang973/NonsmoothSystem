%> test the new coco 
clc
clear
close all
% mu_crit     = 0.025014022367497;
% eta_crit    = 0.095612982101833;
% T_crit      = 4.795364374429504;
% par_SP_Nform        = [-0.1, 0.2, -0.5, 1.781926979010490, 1.6, T_crit, mu_crit, eta_crit];
SN_point_admis  = [-0.1 + 0.2i, -0.1 - 0.2i, -0.5 , 1.781916719621720+0.000010259409,  1.6  , 5.846837671661376]
prob.par  = SN_point_admis;

prob.par2prob  = @par2NForm_Lienard;
% prob.keys   =    keys;
prob.uplim_step = 5000;
prob.ind_p  =    4;
prob.ind_x  =    6;
prob.dir    =    1;
prob.co_dh  =    0.02;
prob.p_span =    [0 1.9];
prob.x_span =    [0 10];
prob.ZeroFunctions = @F_zero;
%> turn on the bifurcation event check 
prob.bifur_detct                = 1;
prob.Events.NumOfEvents         = 6;
prob.Events.MonitorFunctions    = @SN_PD_MonitorFcns;
prob.Events.tag                 = {'SN';'SN';'SN';'PD';'PD';'PD'};
% indexes = [4,  6];
% [t_,Jacob]= get_jacobian(@F_zero, par_deli, indexes)

tic
prob_1        = prob;
prob_2        = prob;
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
plot(u1(4,:), u1(6,:))
hold on
plot(u2(4,:), u2(6,:))
plot(u1(4,SN_locus_1), u1(6,SN_locus_1),'kd','MarkerFaceColor',[1 1 1])
plot(u2(4,SN_locus_2), u2(6,SN_locus_2),'kd','MarkerFaceColor',[1 1 1])
plot(u1(4,PD_locus_1), u1(6,PD_locus_1),'bo','MarkerFaceColor',[1 1 1])
plot(u2(4,PD_locus_2), u2(6,PD_locus_2),'bo','MarkerFaceColor',[1 1 1])
xlabel('r')
ylabel('T')

%> special SN point
SN_points           = u1(:,SN_locus_1)';
% [-0.1 - 0.2i,-0.1 + 0.2i,-0.5,1.78192697900802,1.6,5.84691856641790 ]

SN_point_admis      = SN_points(1,:) ;
[A,B,C,R,T_2_det]   = prob.par2prob(SN_point_admis);
x_00                = IC_generaotr(SN_point_admis(6),R,A,C,-1);
fprintf('Check the SN points! \n')
[Mono_p,Salt_p]     = IC2Floque_Multipliers(SN_point_admis(6),x_00,R,A,C)

%> 
function F = F_zero(par)

        [A,B,C,R, T_2_det]       = par2NForm_Lienard(par);
        F               = T_2_det( par(6) );

end