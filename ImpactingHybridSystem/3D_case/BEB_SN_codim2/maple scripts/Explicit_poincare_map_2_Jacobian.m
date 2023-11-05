%> try to get the numerical return poincare map from analytic formula
Row_11 =  lambda_1^2*(lambda_2 - lambda_3)*exp(T*lambda_1) - lambda_2^2*(lambda_1 - lambda_3)*exp(T*lambda_2) + lambda_3^2*exp(T*lambda_3)*(lambda_1 - lambda_2);
Row_12 =  lambda_1*exp(T*lambda_1)*(lambda_2 - lambda_3) - lambda_2*(lambda_1 - lambda_3)*exp(T*lambda_2) + lambda_3*exp(T*lambda_3)*(lambda_1 - lambda_2);
Row_13 =  (lambda_2 - lambda_3)*exp(T*lambda_1) + (-lambda_1 + lambda_3)*exp(T*lambda_2) + exp(T*lambda_3)*(lambda_1 - lambda_2);

Row_1  = [Row_11, Row_12, Row_13];
%> 
y2 = InitCond(2);
y3 = InitCond(3);
b2 = par_deli(4);
b3 = par_deli(5);
lambda_3 = par_deli(3);
analytic_x_0 = [1,-(b2 *(lambda_1+lambda_2+lambda_3))/(-1+b2)-y2/(-1+b2),-(b3 *(lambda_1+lambda_2+lambda_3))/(-1+b2)-(b3 *y2)/(-1+b2)+y3]
if norm(analytic_x_0' - x_0)/norm(x_0) > 1e-8
    error('There is disagreement between the analytic prediction and the numerical one! ')
end
%> get the partial derivative of the returning condition
d_eq_dT = lambda_1^3*(lambda_2 - lambda_3)*exp(T*lambda_1) - lambda_2^3*...
    (lambda_1 - lambda_3)*exp(T*lambda_2) + lambda_3^3*exp(T*lambda_3)*...
    (lambda_1 - lambda_2) + (lambda_1^2*(lambda_2 - lambda_3)*exp(T*lambda_1)...
    - lambda_2^2*(lambda_1 - lambda_3)*exp(T*lambda_2) + lambda_3^2*exp(T*lambda_3)...
    *(lambda_1 - lambda_2))*y2 +(lambda_1*exp(T*lambda_1)*(lambda_2 - lambda_3) ...
    + (-lambda_1 + lambda_3)*lambda_2*exp(T*lambda_2) + lambda_3*exp(T*lambda_3)...
    *(lambda_1 - lambda_2))*y3;
d_eq_dy2 = lambda_1*exp(T*lambda_1)*(lambda_2 - lambda_3) - lambda_2*(lambda_1 -...
    lambda_3)*exp(T*lambda_2) + lambda_3*exp(T*lambda_3)*(lambda_1 - lambda_2);

d_eq_dy3 = (lambda_2 - lambda_3)*exp(T*lambda_1) + (-lambda_1 + lambda_3)*...
    exp(T*lambda_2) + exp(T*lambda_3)*(lambda_1 - lambda_2);

%> det the dT_dy2 and dT_dy3
dTdy2  = -real(d_eq_dy2/d_eq_dT)
dTdy3  = -real(d_eq_dy3/d_eq_dT)
% dTdy2  = 4.246580303846058; 4.0822
% dTdy3  = 43.151977031641300;  42.9122


%> evaluate the numerator part of the jacobian
Jac_11 =-(b2*lambda_1 + lambda_2 + lambda_3)*lambda_1*(lambda_2 - lambda_3)*(1 + (y2*lambda_1 + lambda_1^2 + y3)*dTdy2)*exp(T*lambda_1) + (lambda_1 - lambda_3)*(b2*lambda_2 + lambda_1 + lambda_3)*(1 + (y2*lambda_2 + lambda_2^2 + y3)*dTdy2)*lambda_2*exp(T*lambda_2) - (b2*lambda_3 + lambda_1 + lambda_2)*(1 + (y2*lambda_3 + lambda_3^2 + y3)*dTdy2)*exp(T*lambda_3)*lambda_3*(lambda_1 - lambda_2);

Jac_12 = -(b2*lambda_1 + lambda_2 + lambda_3)*(1 + lambda_1*(y2*lambda_1 + lambda_1^2 + y3)*dTdy3)*(lambda_2 - lambda_3)*exp(T*lambda_1) + (1 + lambda_2*(y2*lambda_2 + lambda_2^2 + y3)*dTdy3)*(lambda_1 - lambda_3)*(b2*lambda_2 + lambda_1 + lambda_3)*exp(T*lambda_2) - (b2*lambda_3 + lambda_1 + lambda_2)*(1 + lambda_3*(y2*lambda_3 + lambda_3^2 + y3)*dTdy3)*exp(T*lambda_3)*(lambda_1 - lambda_2);

Jac_21 = -((lambda_2 - lambda_3)*lambda_1*(b3*lambda_1 - lambda_2*lambda_3)*exp(T*lambda_1) - (lambda_1 - lambda_3)*lambda_2*(b3*lambda_2 - lambda_1*lambda_3)*exp(T*lambda_2) + exp(T*lambda_3)*(lambda_1 - lambda_2)*lambda_3*(b3*lambda_3 - lambda_1*lambda_2))*(1 + dTdy2);

Jac_22 = -(lambda_2 - lambda_3)*(dTdy3*y2*lambda_1^2 + dTdy3*lambda_1^3 + dTdy3*y3*lambda_1 + 1)*(b3*lambda_1 - lambda_2*lambda_3)*exp(T*lambda_1) + (lambda_1 - lambda_3)*(dTdy3*y2*lambda_2^2 + dTdy3*lambda_2^3 + dTdy3*y3*lambda_2 + 1)*(b3*lambda_2 - lambda_1*lambda_3)*exp(T*lambda_2) - exp(T*lambda_3)*(lambda_1 - lambda_2)*(dTdy3*y2*lambda_3^2 + dTdy3*lambda_3^3 + dTdy3*y3*lambda_3 + 1)*(b3*lambda_3 - lambda_1*lambda_2);

Jac = real([Jac_11, Jac_12; Jac_21, Jac_22]/denom)

eig(Jac)

%> 
fprintf('The indicator of the returning condition! \n')
Row_1*InitCond/denom -1

eig([1.170143904, -1.609566578; 0.831748601, -1.142476193])
