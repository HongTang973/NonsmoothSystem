%> use the CO method to track the solution branch
%> the parameter set in the paper SIADS23
% par                 = [-0.1000 + 0.2000i, -0.1000 - 0.2000i, -0.5 , 1.8,  1.6, 10];

prob.par    =    [-0.1 , 0.2, -0.5 , 1.8,  1.6, 4.791524956512451, 0, 0];

% prob.keys   =    keys;
prob.uplim_step = 2000;
prob.ind_p  =    7;
prob.ind_x  =    [6,8];
prob.dir    =    1;
prob.co_dh  =    1e-2;
prob.p_span = [-0.05 0.05];
prob.x_span = [4 9; -0.6 -0.4];
prob.ZeroFunctions = @F_zero;

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
h1 = plot(u1(3,:)-par_deli(3),u1(4,:)-par_deli(4),'r-','linewidth', 1.5, 'displayname','Addmissible LCO');
hold on
h2 = plot(u2(3,:)-par_deli(3),u2(4,:)-par_deli(4),'k--','linewidth', 1.5,'displayname','Virtual LCO');
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

        [A,B,C,R]       = par2NForm_full(par);
        F               = [Det(A,R,C,par(6) );
            root2tangency(A,R,C,  par(6) )];

end

%

%> zero function condition
function f = Det(A,R,C,T)
f = det(R*expm(A*T) - eye(length(C)));

end

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

%> derivative for single variable function
function varargout = delicate_root(tspan, par,index, var)
par(index)  = var;
[A,B,C,R]   = par2NForm(par);
fs          = 1e6;
delta       = 1/fs;
T           = tspan(1):delta:tspan(2);
%
Det         = @(t) det(R*expm(A*t) - eye(length(C)));

%> initialize the container of roots
% roots_ind   = [];
%> set the default sign as positive

fval_list   = zeros(1, length(T));
%> calculate the value of the det function and detect the sign change
for i=1:length(T)
    fval_list(i) = Det(T(i));
end

[fval, ind]  = max(fval_list);
%> output the results at the steady status with given valve lift
if nargout >= 0; varargout{1} =  fval;       end
if nargout >= 2; varargout{2} =  T(ind);            end
end