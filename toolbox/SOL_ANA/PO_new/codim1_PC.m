
% function [u,iter] = codim1_PC(ZeroFunctions,par,ind_p,ind_x,p_span,x_span,dir)

% Parameter Analysis of the found period one orbit

% Self coded by knowledg from the book
% Introduction to Numerical Continuation Methods --- Allgower

% ind_p:            specify the varying parameter
% ind_x:            specify the solution state 
% ZeroFunctions:    function handle of zero funcitons of OP
% par:              give the para meter vector
% p_span:           specify the range of parameter to be continued
% dir:              specify the direction of the solution branch
function [u,iter] = codim1_PC(prob)
% 
ind_p               = prob.ind_p;
ind_x               = prob.ind_x;
ZeroFunctions       = prob.ZeroFunctions;
par                 = prob.par;
p_span              = prob.p_span;
x_span              = prob.x_span;
dir                 = prob.dir;
% now par must be a vector now
% arc length stepsize
default_h = 0.001;

dh = min(default_h,(p_span(2)-p_span(1))/200);

% DETERMINE THE LOCATION of the parameter's partial derivative 
index       = [ind_p,ind_x]; 
prob.index  = index;

if  size(x_span,1)~= length(ind_x) || size(p_span,1)~= length(ind_p)
    warndlg('Wrong length match of p and x.')
end

% searching direction of parameter space
% dir = dir; 


%   **************** Codimension one ******************************    %


% step 1: set the closed form zero equations
% H(x,a) = 0
% u_0 = [x,a];
% H(u_0) = 0;

% step 2:  using the local manifold information to span the solution branch
%  Jacobian matrix : H'(u_0)
%  Rank(H'(u_0)) = n;  the one-dimensional kernel ker(H')
%  if u(a) denotes the solution branch
%               H'(u) u'(a) = 0
%  use the condition of det([H'; u']) > 0 to choose the right tangent
%  vector


% step 3: stepsize of prediction by the arclength
% ds =  norm(du_da) da



% step 4: Euler predictor


% step 5: corrector step

%  ---------------------       ----------------------------      -------- %

% ##### Step 1: define zero functions

  % scalar function 1D searching



% check the par is a vector
if size(par,2)==1;  else;  par = par.';  end

% initial condition of the continuation system

u = par; % container of the  parameter track

iter = [];



[prob_det] = ZeroFunctions (prob);

if norm(prob_det) > 5*1e-6
        warndlg( 'Not zero parameter point!' )
end
% 

% ##### Step 2: get the Jacobian 
% par = {'U', 'r', 'damp_1','damp_2', 'damp_3', 'T'}



[t_,~] = get_jacobian(prob); 

trans_flag =  norm(t_);

reg = u(ind_p,end)>p_span(1) && u(ind_p,end)<p_span(2)...
     && all(u(ind_x,end)-x_span(:,1)>=0) && all(x_span(:,2)-u(ind_x,end)>=0);





new_par = u(:,end); dis2init = 1;
while trans_flag && reg
    
   % first predictor
    % using the tangent vector produced by the Jacobian
    new_par = Predictor(new_par,index,ind_p,t_,dh,dir,dis2init);
    
    prob_temp               = prob;
    prob_temp.par           = new_par;
 
    % corrector: now the newly predicted point in parameter space, but this should be corrected by iteration
    num = 1;
    [u,iter,dis2init,iscirc] = corrector(prob,u,num,iter);
    
    prob_temp.par            =  u(:,end);
    [t_,~]                   = get_jacobian(prob_temp); 
    
    trans_flag =  norm(t_);
    
 
  reg = iscirc && u(ind_p,end)>p_span(1) && u(ind_p,end)<p_span(2)...
     && all(u(ind_x,end)-x_span(:,1)>=0) && all(x_span(:,2)-u(ind_x,end)>=0);
    
end

function new_par = Predictor(new_par,index,ind_p,t_,dh,dir,dis2init)
    
    if abs(t_(1))<1e-2  % turning point 
        % the index =1, since we put the parameter in the first place when
        % we get our jacobian
        
        dh = min(0.02*dh, 0.1*abs(dis2init(1)));
        K = null(t_');
        co = [1;rand(length(t_)-1,1)];
        co = co/norm(co);
        % new searching direction 
        t_ = [t_,K]*co;
           
    elseif  abs(dis2init(1))<2*dh && sign(dis2init(1)*t_(1))<0 && dis2init(2)<1e-3 
        % the second condition is the loop returning conditon to prevent the
        % low starting speed at first steps
        dh = min(dh*0.1,0.5*abs(dis2init(1))*new_par(ind_p));
        
    else
        % plain predictor step
        
    end
    new_par(index) = new_par(index) + dir * dh* t_;
 %%   
 function [u,iter,dis2init,iscirc] = corrector(prob,u,num,iter)
     ind_p   = prob.ind_p;
     ind_x   = prob.ind_x;
     index   = [ind_p, ind_x]; 
     par_0   = u(:,1); % the first point 
     ZeroFunctions = prob.ZeroFunctions;
     prob_temp     = prob;
     new_par       = prob_temp.par;
     while 1
         [prob_det]      = ZeroFunctions(prob_temp);
         
         [~,J]           = get_jacobian(prob_temp);
         
         alt             = -J'*inv(J*J')*prob_det;
         
         new_par(index)  = new_par(index) + alt;
         
         prob_temp.par   = new_par;
         
         num             = num + 1;
         
         
         % check convergence
         if norm(prob_det) <= 1e-12 || num >=100
             
             u          = [u,new_par];
             iter       = [iter,num];
             
             
             % incase a closed circle solution
             dis2init   = [(new_par(ind_p)-par_0(ind_p))/norm(par_0(ind_p));
                 norm(new_par(ind_x)-par_0(ind_x))/norm(par_0(ind_x))] ;
             
             iscirc     = norm(new_par(index)-par_0(index))/norm(par_0)>1e-6;
             
             if num >=100
                 disp('Iteration times limited!!')
             end
             
             
             %
             break;
         end
     end