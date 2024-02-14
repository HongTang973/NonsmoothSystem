%> ------------- built in funcitons -------------------------- %

function new_par = Predictor(prob,last_par,t_,dh,dis2init)
ind_p               = prob.ind_p;
ind_x               = prob.ind_x;
index               = [ind_p,ind_x];
dir                 = prob.dir;
if abs(t_(1))<1e-2  % turning point
    % the index =1, since we put the parameter in the first place when
    % we get our jacobian
    fprintf('The turning point is detected!\n')
    dh = max(prob.min_h, min(0.02*dh, 0.1*abs(dis2init(1))));
    K = null(t_');
    co = [1;rand(length(t_)-1,1)];
    %> strategy 1
    co = co/norm(co);
    % new searching direction
    t_ = [t_,K]*co;
    
elseif  abs(dis2init(1))<2*dh && sign(dis2init(1)*t_(1))<0 && dis2init(2)<1e-3
    % the second condition is the loop returning conditon to prevent the
    % low starting speed at first steps
    % dh = min(dh*0.1,0.5*abs(dis2init(1))*new_par(ind_p));
    dh = min(dh, 0.001);
    
else
    % plain predictor step
    
end
new_par         = last_par;
new_par(index)  = last_par(index) + dir * dh* t_/norm(t_);