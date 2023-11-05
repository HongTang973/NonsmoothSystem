function  [t_,Jacob]= get_jacobian(ZeroFunctions, par, index)
% ZeroFunctions: define the problem 
% par: the specific par point 

if isvector(par)
    par     =    par;
elseif isstruct(par)
    par     =    cell2mat(struct2cell(par));
else
    warndlg('Wrong input in the codim1_PC function!')
end
% 
value = ZeroFunctions(par);

% if norm(value)>1e-2
%     fprintf('It is too far from the equilibrium')
% end

    
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

for i = 1:n_vp

    par1 = par;
    par2 = par;
    delta = init_delta;
    tol = 1;
    
    par1(index(i)) = par(index(i)) - delta;
    par2(index(i)) = par(index(i)) + delta;
   
    % initial step
    Jac_temp = (ZeroFunctions (par2)-ZeroFunctions (par1))/2/delta;
    
    while  delta > min_delta
    %  iteration step
    delta   = delta /2;
    par1(index(i)) = par(index(i)) - delta;
    par2(index(i)) = par(index(i)) + delta;

    Jac_new = (ZeroFunctions (par2)-ZeroFunctions (par1))/2/delta;
    
    tol = norm(Jac_new - Jac_temp)/norm(Jac_temp);
    
    Jac_temp = Jac_new;
     
    if tol < 1e-6
        Jacob(:,i) = 0.5*(Jac_new + Jac_temp);
        break;
    elseif tol<1e-5
        Jacob(:,i) = 0.5*(Jac_new + Jac_temp);
    elseif tol<1e-4
        Jacob(:,i) = 0.5*(Jac_new + Jac_temp);
    elseif tol<1e-2
        Jacob(:,i) = 0.5*(Jac_new + Jac_temp);
    end
    
    end 
end
if size(Jacob,1) ~= size(Jacob,2)
    % get the tangent vector induced by the Jacobian
    if any(isnan(Jacob))
        keyboard
    else
        t_ = null(Jacob);
    end
    %>
    try
        det_dir = det([Jacob;t_']);
    catch ME
        keyboard
    end
    %>
    if sign(det_dir)==1
        
    elseif sign(det([Jacob;t_']))==-1
        t_ = -t_;
    else
        t_ = 0;
    end
else
    t_ = 0;
end

end