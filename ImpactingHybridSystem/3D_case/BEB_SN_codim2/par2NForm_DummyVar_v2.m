function [A,B,C,R,T_2_det] = par2NForm_DummyVar_v2(par)
% this version of adding mu dependence is after the comments form David
% Simpson
% F = A*x + M*mu + A_1*mu*x

% par = {'lambda_1', 'lambda_2','lambda_3', 'b2','b3', 'T', 'mu', 'eta'}
%> append the system par with two dummy variables to do the analysis

% par =  [-0.1, 0.2, -0.5, 1.781926979010490, 1.6, 5.846868671661380, 0, 0];

%> use the normal form parameters to form the matrixes
% par(1) = par(1) + par(8);
%par(3) = par(3) - par(7);
par(4)  = par(4) + par(8);
% par(5) = par(5) + par(7);

A_1 = []
[A,B,C,R,T_2_det] = par2NForm_Lienard(par);
A = A + par(7)*A_1;
end