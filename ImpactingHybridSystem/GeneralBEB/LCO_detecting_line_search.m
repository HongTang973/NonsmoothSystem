%>  ================================================================ 
%> @brief  Use structure type to construnct the OOP valve+pipe model for gas
%> 
%> Initialize the numerical tool
%>  
%> @param p input structure type which could be empty
%>
%> @retval p the constructed structure type to store the information of the Valve  model and the functions
%  =================================================================
function [sign,LOCI,Max,vector,F_1,COND] = LCO_detecting_line_search(T,R,EA,sign_V,C,equi_type)
 
Mono_A = EA(T);
M = R*Mono_A- eye(length(C));
f_1 = det(M);
PA = R*Mono_A;
% construct a transform matrix
Trans = eye(length(C));
Trans(C>0,:) = []; % get a (n-1) x n tranform matrix
%
ROW_M  = C*M; % to check the vadality of the remaining equation
COL_M = Trans*M(:,C>0);
%
M(C>0,:) =[];
M(:,C>0) =[];
% detect the singular point
COND = rcond(M);
Ans = -M\COL_M*equi_type;
% substitute the sub-system's solution to the remaining eqaution
%
try
F_1 = ROW_M*equi_type*(det(M)*C' - Trans'*adjugate(M)*COL_M);
catch msg
 F_1 = ROW_M*(equi_type*C'+ Trans'*Ans);
 F_1 = f_1;
end
% F_1 = ROW_M*(equi_type*C'+ Trans'*Ans);
% calculate the eigvalue of the composed matrix first get the tiny element
% vanished
PA(abs(PA)<1e-8)=0;
[JV,JD]=eig(PA);
% select just real eigenvalue
eig_A = diag(JD);
eig_A = eig_A(eig_A==real(eig_A));
JV = JV(:,eig_A==real(eig_A));
% selelct real eigenvector
JV = JV(:,sum(JV==real(JV))==length(C));
eig_A = eig_A(sum(JV==real(JV))==length(C));
%
[MaxEig,~] = max(abs(eig_A));
[~,ind2] = min(abs(abs(eig_A)-1));
MaxEig = max(MaxEig,max(abs(diag(JD))));
V=JV(:,ind2);
if  isempty(V)
    Max = max(abs(diag(JD)))-1;
    LOCI= sort((diag(JD)))-1;
    sign =-1;
    vector = zeros(length(C),1);
else
    
    if V(C>0)
        V = V/V(C>0);
    end
    vector =equi_type*V;
    Max = MaxEig-1;
    LOCI= sort((diag(JD)))-1;
    sign =sign_V(vector);
end
%

function B = adjugate(A)
% This finds the adjugate (adjoint) of square matrix A,
% and is valid even if A is singular or complex-valued.
% With u, s, and v obtained from [u,s,v] = svd(A), it
% makes use of the identity adj(A) = det(u*v')*v*adj(s)*u',
% which holds even if A and s are singular.  The expression,
% diag(prod(reshape(s0(ix),n-1,n),1)), accomplishes the
% evaluation of adj(s), each of whose diagonal elements
% is the product of all but one of the diagonal elements
% of s.  This requires that A be n x n where n >= 2.
% Roger Stafford - 10/18/06
[m,n] = size(A);
if (m ~= n) | (n < 2)
 error('Matrix A should be size n x n with n >= 2.')
end
[u,s,v] = svd(A);
s0 = diag(s);
ix = toeplitz(ones(n-1,1),[1 zeros(1,n-1)]) ...
     + repmat((1:n-1)',1,n);
B = det(u*v')*v*diag(prod(reshape(s0(ix),n-1,n),1))*u';
% ---
end
end