function J = Shilnikov_IHS_DFDX(x, p, mode) %#ok<INUSD>
global const
%> predefined: the term of mu^2*A00*x
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

%>  'hspo'-compatible encoding of the Jacobian w.r.t. state variables.
%>  par = {'rho', 'omega', 'lambda','r0', 'c0', 'mu','eta'}
rho     = p(1);
omega   = p(2);
lambda  = p(3);
mu      = p(6);
eta     = p(7);
%
x1 = x(1,:);
x2 = x(2,:);
x3 = x(3,:);
%
J = zeros(3,3,numel(x1));
J(1,1,:) = a11*mu^2     + 2*c11*x1      + c12*x2    + c31*x3    + rho;
J(1,2,:) = a12*mu^2     + c12*x1        + 2*c22*x2  + c23*x3    + omega;
J(1,3,:) = a13*mu^2     + c31*x1        + c23*x2    + 2*c33*x3;
%
J(2,1,:) = a21*mu^2     + 2*d11*x1      + d12*x2    + d31*x3    - omega;
J(2,2,:) = a22*mu^2     + d12*x1        + 2*d22*x2  + d23*x3    + rho;
J(2,3,:) = a23*mu^2     + d31*x1        + d23*x2    + 2*d33*x3  + 1;
%
J(3,1,:) = a31*mu^2     + 2*e11*x1      + e12*x2    + e31*x3;
J(3,2,:) = a32*mu^2     + e12*x1        + 2*e22*x2  + e23*x3;
J(3,3,:) = a33*mu^2     + e31*x1        + e23*x2    + 2*e33*x3  + lambda;
%
end
