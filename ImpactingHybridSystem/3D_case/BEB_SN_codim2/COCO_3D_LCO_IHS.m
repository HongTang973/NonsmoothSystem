function [data, prob_det] = COCO_3D_LCO_IHS (prob, data, par)

% par = {'lambda_1', 'lambda_2','lambda_3', 'b2','b3', 'T'}
% par =  [-1, 0.2, -1, 2.5, 0.5, 0.875926384431162];
%%  define the linear part & get the matrix A  
par2prob  = prob.par2prob;
[A,B,C,R] = par2prob(par);
% sign_V = @(Y) sign(C*A_1*Y);
% %
% A_s = (eye(8)-(B*C*A)/(C*A*B))*A;
T = par(6);
%
Mono_A = expm(A*T);
% Mono_A = EA(T);
M = R*Mono_A - eye(length(C));
%
prob_det = det(M);
% penalty function to escape the zero value line
% if  abs(lambda_1)<1e-4 
%     
%     prob_det = prob_det + 1/lambda_1;
%     
% elseif abs(lambda_2)<1e-4 
%     
%     prob_det = prob_det + 1/lambda_2;
%     
% elseif abs(lambda_3)<1e-4 
%     
%     prob_det = prob_det + 1/lambda_3;
%     
% else
%     % 
% end


end