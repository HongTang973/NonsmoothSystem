function [A,B,C,R,T_2_det] = par2NForm_DummyVar_David(par)
% par = {'lambda_1', 'lambda_2','lambda_3', 'b2','b3', 'T', 'mu', 'eta'}
%> append the system par with two dummy variables as to track the codim2 curve 
%> to do the analysis

%> e.g. this is a codimension 2 point for BEB-SN
% par =  [-0.1, 0.2, -0.5, 1.781926979010490, 1.6, 5.846868671661380, 0, 0];

%> use the normal form parameters to form the matrixes
% par(1) = par(1) + par(8);
% par(3) = par(3) - par(7);

par(4) = par(4) + par(8);
% par(5) = par(5) + par(7);

% replace the mu dependence as david wanted, which will be mu scaled
% perturbation of the Jacobian matrix
[A,B,C,R,T_2_det] = par2NForm_Lienard_David_CD2(par);
end