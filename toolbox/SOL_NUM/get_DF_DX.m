
%>  ================================================================ 
%> @brief  General method to numerically evaluate the jacobian matrix
%> 
%> 
%>  
%> @param F the function handle of the customed function
%> @param y the point to around which the jacobian will be evaluated
%> @param index specify the variables to vary to get the partial derivatives to make up the Jacobian
%> @param p the structure or the class of the physical object <if needed>
%>
%> @retval Jacob numerical Jacobian
%  =================================================================
function  [Jacob]= get_DF_DX(varargin)
if nargin  == 3
    F       = varargin{1};
    y       = varargin{2};
    index   = varargin{3};
    p       = [];
    func    = @(t,y,prob) F(t,y);
elseif nargin == 4
    F       = varargin{1};
    y       = varargin{2};
    index   = varargin{3};
    p       = varargin{4};
    func    = @(t,y,prob) F(t,y,prob);
else
    error('Wong input!')
end
t =0;
value = func(t,y,p);

if ~isreal(value)
    keyboard
end
% define the number of equations
n_eq        = length(value);
% define the number of varying parameters
n_vp        = length(index);
% allocate memory for Jacobian
Jacob       = zeros(n_eq,n_vp);
init_dh     = 1e-3;
min_dh      = 1e-12;
% h  = 1.0e-8*( 1.0 + abs(par(index(end))) );
% hi = 0.5./h;

for i = 1:n_vp
    y1 = y;
    y2 = y;
    dh = init_dh;
    y1(index(i)) =y(index(i)) - dh;
    y2(index(i)) =y(index(i)) + dh;
    Jac_temp     = (func(t,y2,p) - func(t,y1,p))/2/dh;
    while  dh > min_dh
        dh       = dh /10;
        y1(index(i)) =y(index(i))   - dh;
        y2(index(i)) =y(index(i))   + dh;
        Jac_new      = (func(t,y2,p)   - func(t,y1,p))/2/dh;
        tol          = norm(Jac_new - Jac_temp)/norm(Jac_temp);
        Jacob(:,i)   = 0.5*(Jac_new + Jac_temp);
        Jac_temp     = Jac_new;
        % terminate if the tollerance is arrived
        if tol < 1e-9
            break;
        end
    end
end


end