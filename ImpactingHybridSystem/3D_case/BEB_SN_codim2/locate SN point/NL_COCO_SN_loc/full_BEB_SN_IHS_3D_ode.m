%> from parameter to the shilnikov case autonomous system
function dxdt = full_BEB_SN_IHS_3D_ode(x, par, mode)
% 
global const
% 
[A,B,C,R]       = par2NForm_DummyVar(par);
mu              = par(7);
eta             = par(8);
M               = [0;0;-1];
% 
A00             = const.A00;
B00             = const.B00;
C00             = const.C00;
% 
% in matrix form 
% NL_term = mu^2*A00*[x1;x2;x3] + B00*[x1*x2; x2*x3; x3*x1] + C00*[x1^2; x2^2; x3^2]
% --- index = [2 3 1]
% -- bilinear term:     y.*y(index)
% -- quadratic term:    y.*y
index                = [2 3 1];
dxdt (:,:)           = A*x + M*mu + 0*mu^2*A00*x + B00*(x.*x(index,:)) + C00*x.^2;
end