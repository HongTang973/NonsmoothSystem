%%
function [u,iter,dis2init,iscirc,discre_err] = corrector(prob,u, new_par,num,iter)
ind_p               = prob.ind_p;
ind_x               = prob.ind_x;
ZeroFunctions       = prob.ZeroFunctions;
tol     = 1e-12;
index   = [ind_p,ind_x];
par     = u(:,1); % the first point
while 1
    [prob_det]  = ZeroFunctions (new_par);

    [~,J]       = get_jacobian(ZeroFunctions, new_par, index);

    alt         = -real(J'*inv(J*J')*prob_det);
    if norm(alt)>1e1 %> avoid singular prediction
        % keyboard
        alt     = alt./norm(alt).*prob.co_dh;
        tol     = 1e-6;
    end
    new_par(index) = new_par(index) + alt;

    num         = num + 1;


    % check convergence
    discre_err = norm(prob_det);
    if  discre_err <= tol || num >=100

        u = [u,new_par];
        iter = [iter,num];
        % incase a closed circle solution
        dis2init = [(new_par(ind_p)-par(ind_p))/(1+norm(par(ind_p)));
            norm(new_par(ind_x)-par(ind_x))/norm(par(ind_x))] ;

        iscirc = norm(new_par(index)-par(index))/(1+norm(par(index)))>1e-12;

        if num >=100
            fprintf('Iteration times limited!! with %d\n',norm(prob_det))
        end

        %
        break;
    end
end