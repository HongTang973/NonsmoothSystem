% function  [t_,Jacob]= get_jacobian(ZeroFunctions, par, index)
function  [t_,Jacob]= get_jacobian(prob)
    %  prob: the problem of solution 

% ZeroFunctions: define the problem 
% par: the specific par point 
ZeroFunctions = prob.ZeroFunctions;
par           = prob.par;
index         = prob.index;
% 
value = ZeroFunctions(prob);

    
% define the number of equations     
n_eq = length(value);
% define the number of varying parameters
n_vp = length(index);
% allocate memory for Jacobian
Jacob = zeros(n_eq,n_vp);

init_delta = 1e-3; 
min_delta = 1e-12;


% h  = 1.0e-8*( 1.0 + abs(par(index(end))) );
% hi = 0.5./h;
prob_tmp_1 = prob;
prob_tmp_2 = prob;

for i = 1:n_vp

    par1 = par;
    par2 = par;
    delta = init_delta;
    tol = 1;
    
    par1(index(i)) = par(index(i)) - delta;
    par2(index(i)) = par(index(i)) + delta;
    
    % initial step
    prob_tmp_1.par = par1;
    prob_tmp_2.par = par2;

    Jac_temp = (ZeroFunctions (prob_tmp_1)-ZeroFunctions (prob_tmp_2))/2/delta;
    
    while  delta > min_delta
    %  iteration step
    delta = delta /10;
    par1(index(i)) = par(index(i)) - delta;
    par2(index(i)) = par(index(i)) + delta;

    prob_tmp_1.par = par1;
    prob_tmp_2.par = par2;

    Jac_new = (ZeroFunctions (prob_tmp_1)-ZeroFunctions (prob_tmp_2))/2/delta;
    
    tol = norm(Jac_new - Jac_temp)/norm(Jac_temp);
    
    Jacob(:,i) = 0.5*(Jac_new + Jac_temp);
    
    Jac_temp = Jac_new;
    
    if tol < 1e-6
        break;
    end
    
    end 
end

% get the tangent vector induced by the Jacobian 
t_ = null(Jacob);

if sign(det([Jacob;t_']))==1
    
elseif sign(det([Jacob;t_']))==-1
    t_ = -t_;
else
    t_ = 0;
end

end