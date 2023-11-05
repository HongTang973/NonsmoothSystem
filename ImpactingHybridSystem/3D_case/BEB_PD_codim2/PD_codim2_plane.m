%> use the CO method to track the solution branch
%> the parameter set in the draft paper of Shilnikov 
% par                 = [0.1, 1, 0.3, 0.640514496669177, 0.788, 9.95899559243775];
clc
prob.par        =    [0.1, 1, 0.3, 0.640514496669177,0.788,9.95899594243774];
par_deli        =    [0.1, 1, 0.3, 0.640514496669177,0.788,9.95899594243774];
% prob.keys   =    keys;
prob.uplim_step = 2000;
prob.ind_p  =    4;
prob.ind_x  =    [5,6];
prob.dir    =    1;
prob.co_dh  =    1e-2;
prob.p_span =    [0.9*par_deli(4) 1.05*par_deli(4)];
prob.x_span =    [0 1; 5 12];
prob.ZeroFunctions = @F_zero;
indexes = [4, 5, 6];
[t_,Jacob]= get_jacobian(@F_zero, par_deli, indexes)
%> 
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
x_t_range_1 = [-0.05:0.001:0.05];
x_t_range_2 = [-0.02:0.0001:0];
FIG1 = figure;
h1 = plot(u1(4,:)-par_deli(4),u1(5,:)-par_deli(5),'r-','linewidth', 1.5, 'displayname','Addmissible LCO');
hold on
h2 = plot(u2(4,:)-par_deli(4),u2(5,:)-par_deli(5),'k--','linewidth', 1.5,'displayname','Virtual LCO');
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


%> 
function F = F_zero(par)

        [A,B,C,R,T_2_det]       = par2NForm_Shilnikov(par);
        F                       = [T_2_det( par(6) );
            root2tangency(A,R,C,  par(6) )];

end

%



%> the tangency condition
function [ diff_det] = root2tangency(A,R,C,T)
Det         = @(t) det(R*expm(A*t) - eye(length(C)));
init_delta  = 1e-3;
min_delta   = 1e-12;
dt          = init_delta;
d_tmp           = 10;
iter            = 0;
while dt > min_delta
    d_new      = ( Det(T+dt) - Det(T-dt))/2/dt;
    tol        = abs(d_new - d_tmp);
    dt         = dt /2;
    d_tmp      = d_new;
    iter       = iter + 1;
    if tol < 1e-9
        break;
    end
    
end
diff_det = d_new;
end

