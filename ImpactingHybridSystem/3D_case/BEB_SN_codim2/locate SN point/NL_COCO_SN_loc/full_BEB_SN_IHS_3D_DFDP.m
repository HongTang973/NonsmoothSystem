function J = Shilnikov_IHS_DFDP(x, p, mode)
global const
%> predefined: the term of mu^2*x
a11 = const.A00(1,1);
a12 = const.A00(1,2);
a13 = const.A00(1,3);
a21 = const.A00(2,1);
a22 = const.A00(2,2);
a23 = const.A00(2,3);
a31 = const.A00(3,1);
a32 = const.A00(3,2);
a33 = const.A00(3,3);
%> predefined  bilinear terms
c12 = const.B00(1,1);
c23 = const.B00(1,2);
c31 = const.B00(1,3);
%
d12 = const.B00(2,1);
d23 = const.B00(2,2);
d31 = const.B00(2,3);
%
e12 = const.B00(3,1);
e23 = const.B00(3,2);
e31 = const.B00(3,3);
%> predefined  quadratic terms
c11 = const.C00(1,1);
c22 = const.C00(1,2);
c33 = const.C00(1,3);
%
d11 = const.C00(2,1);
d22 = const.C00(2,2);
d33 = const.C00(2,3);
%
e11 = const.C00(3,1);
e22 = const.C00(3,2);
e33 = const.C00(3,3);

%

%   'hspo'-compatible encoding of the Jacobian w.r.t. state variables.
%> par = {'rho', 'omega', 'lambda','r0', 'c0', 'mu','eta'}
rho     = p(1);
omega   = p(2);
lambda  = p(3);
mu      = p(6);
eta     = p(7);
%
x1      = x(1,:);
x2      = x(2,:);
x3      = x(3,:);
%
J = zeros(3,8,numel(x1));
%
J(1,1,:) = x1;
J(1,2,:) = x2;
J(1,7,:) = 2*a11*mu*x1 + 2*a12*mu*x2 + 2*a13*mu*x3;
%
J(2,1,:) = x2;
J(2,2,:) = -x1;
J(2,7,:) = 2*a21*mu*x1 + 2*a22*mu*x2 + 2*a23*mu*x3;

%
J(3,3,:) = x3;
J(3,7,:) = 2*a31*mu*x1 + 2*a32*mu*x2 + 2*a33*mu*x3 + 1;
%


end
