
function [u,iter] = codim1_PC_skip0(ZeroFunctions,par,ind_p,ind_x,p_span,x_span,dir)

% Parameter Analysis of the found period one orbit

% Self coded by knowledg from the book
% Introduction to Numerical Continuation Methods --- Allgower

% ind_p:            specify the varying parameter
% ind_x:            specify the solution state 
% ZeroFunctions:    function handle of zero funcitons of OP
% par:              give the para meter vector
% p_span:           specify the range of parameter to be continued
% dir:              specify the direction of the solution branch


% arc length stepsize

dh = min(0.001,(p_span(2)-p_span(1))/200);

% DETERMINE THE LOCATION of the parameter's partial derivative 
index = [ind_p,ind_x]; 

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



[prob_det] = ZeroFunctions (par);

if norm(prob_det) > 5*1e-6
        warndlg( 'Not zero parameter point!' )
end
% 

% ##### Step 2: get the Jacobian 
% par = {'U', 'r', 'damp_1','damp_2', 'damp_3', 'T'}



[t_,~] = get_jacobian(ZeroFunctions, par, index); 

trans_flag =  norm(t_);

% the flag to continue or stop the continuation
reg = u(ind_p,end)>p_span(1) && u(ind_p,end)<p_span(2)...
     && all(u(ind_x,end)-x_span(:,1)>=0) && all(x_span(:,2)-u(ind_x,end)>=0);





new_par = u(:,end); dis2init = 1; n_skip =1;
while trans_flag && reg
    
    % first predictor
    % using the tangent vector produced by the Jacobian
    new_par = Predictor(new_par,index,ind_p,t_,dh,dir,dis2init);
    
    
    % corrector
    num = 1;
    [u,iter,dis2init,iscirc] = corrector(ZeroFunctions,u, new_par, ind_p, ind_x,num,iter);
        
       
%         [prob_det] = ZeroFunctions (new_par);
%         
%         [~,J] = get_jacobian(ZeroFunctions, new_par, index);
%         
%         alt = -J'*inv(J*J')*prob_det;
%         
%         new_par(index) = new_par(index) + alt;
%         
%         num = num + 1;
%         
%         
%         % check convergence
%         if norm(prob_det) <= 1e-9 || num >=100
%             
%             u = [u,new_par];
%             iter = [iter,num];
%             
%             
%             % incase a closed circle solution
%             dis2init = [(new_par(ind_p)-par(ind_p))/norm(par(ind_p));
%                norm(new_par(ind_x)-par(ind_x))/norm(par(ind_x))] ;
%            
%             iscirc = norm(new_par(index)-par(index))/norm(par)>1e-6;
%                
%             if num >=100
%                 disp('Iteration times limited!!')
%             end
%             
%             
%             % 
%             break;
%         end
        
        % 
        
        
        
        %   modified by peter on 18/8/2022: when close to the zero parameter
        %     value, skip the point with another new starting  point
        if any(abs(new_par(index))<1e-3)
            % skip the zero value to give a new initial point
            t_temp = new_par(index)- par(index);
            new_par(index) = new_par(index) + 0.1*t_temp/norm(t_temp);
            
%             new_par = Predictor(new_par,index,ind_p,t_,dh,dir,dis2init);
%             
%             ZeroFunctions(new_par)
            [u_temp,~,~,iscirc] = corrector(ZeroFunctions,u, new_par, ind_p, ind_x,num,iter);
            
            u = [u,u_temp(:,end)];
            iter = [iter,-n_skip]; % 0 as a sign of skip
            n_skip = n_skip+1;
        end
        
            
    %  judge the direction of the searching direction of current point
    [t_,~] = get_jacobian(ZeroFunctions, u(:,end), index); 
    
    trans_flag =  norm(t_);
    
 
    
  reg = iscirc && u(ind_p,end)>p_span(1) && u(ind_p,end)<p_span(2)...
     && all(u(ind_x,end)-x_span(:,1)>=0) && all(x_span(:,2)-u(ind_x,end)>=0);
    
end
%% 
function new_par = Predictor(new_par,index,ind_p,t_,dh,dir,dis2init)
    
    if abs(t_(1))<1e-2  % turning point 
        % the index =1, since we put the parameter in the first place when
        % we get our jacobian
        
        dh = min(dh, 0.1*abs(t_(ind_p)));
        K = null(t_');
        co = [1;rand(length(t_)-1,1)];
        co = co/norm(co);
        % new searching direction 
        t_ = [t_,K]*co;
           
    elseif  abs(dis2init(1))<2*dh && sign(dis2init(1)*t_(1))<0 && dis2init(2)<1e-3 
        % the second condition is the loop returning conditon to prevent the
        % low starting speed at first steps
        dh = min(dh*0.01,0.5*abs(dis2init(1))*new_par(ind_p));
        
    else
        % plain predictor step
        
    end
    new_par(index) = new_par(index) + dir * dh* t_;
    
  %%   
 function [u,iter,dis2init,iscirc] = corrector(ZeroFunctions,u, new_par, ind_p, ind_x,num,iter)
     index = [ind_p,ind_x]; 
     par =u(:,1); % the first point 
     while 1
         [prob_det] = ZeroFunctions (new_par);
         
         [~,J] = get_jacobian(ZeroFunctions, new_par, index);
         
         alt = -J'*inv(J*J')*prob_det;
         
         new_par(index) = new_par(index) + alt;
         
         num = num + 1;
         
         
         % check convergence
         if norm(prob_det) <= 1e-9 || num >=100
             
             u = [u,new_par];
             iter = [iter,num];
             
             
             % incase a closed circle solution
             dis2init = [(new_par(ind_p)-par(ind_p))/norm(par(ind_p));
                 norm(new_par(ind_x)-par(ind_x))/norm(par(ind_x))] ;
             
             iscirc = norm(new_par(index)-par(index))/norm(par)>1e-6;
             
             if num >=100
                 disp('Iteration times limited!!')
             end
             
             
             %
             break;
         end
     end