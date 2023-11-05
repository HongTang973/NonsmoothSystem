%>  ================================================================ 
%> @brief  return the det value to assert the critical condition of LCO's
%> exixtence
%> 
%> Form the det problem from the given parameters of the system
%>  
%> @param  T the guessed return time
%> @param  R the reset map
%> @param  A the Jacobian of the vector flow field near the BEB
%> @param  C the observation vector
%> @param  equi_type indicator of the equlibrium type, either admissible or
%> virtual
%>
%> @retval p the constructed structure type to store the information of the Valve  model and the functions
%  =================================================================
function [sign_V,LOCI,Max,vector,F_1,COND] = LCO_Det_search(T,R,A,C,equi_type)
method  = true;
%> monodromy matrix construction
Mono_A = expm(A*T);
%> the determinant condition
M = R*Mono_A- eye(length(C));
%> the straight way to calculate the determinant of the characteristic
%> matrix
f_1 = det(M);

PA = R*Mono_A;

% construct a transform matrix in case  random nonzero component in C vector
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
if method
    F_1 = f_1;
else
   
    Ans = -M\COL_M*equi_type;
    % substitute the sub-system's solution to the remaining equation
    %> formulate the candidate of initial condition 
    IC  = (equi_type*C'+ Trans'*Ans);
    try
        F_1 = ROW_M*equi_type*(det(M)*C' - Trans'*adjugate(M)*COL_M);
    catch msg
        F_1 = ROW_M*IC;
    end
    
end

% calculate the eigvalue of the composed matrix first get the tiny element
% vanished
PA(abs(PA)<1e-9)=0;
%> solve the eigen values and the eigen vectors of the map
if any(any(isnan(PA))) || any(any(~isfinite(PA)))
    keyboard
end
[JV,JD]=eig(PA);

% select just real eigenvalue and corresponding eigenvectors
eig_A = diag(JD);
eig_A = eig_A(abs(imag(eig_A))<1e-9);
%> 
JV = JV(:,eig_A==real(eig_A));
% selelct real eigenvector
JV = JV(:,sum(JV==real(JV))==length(C));
eig_A = eig_A(sum(JV==real(JV))==length(C));
%> get the eigenvalue with the greatest magnitude
[MaxEig,~] = max(abs(eig_A));
[~,ind2] = min(abs(abs(eig_A)-1));
MaxEig = max(MaxEig,max(abs(diag(JD))));
V=JV(:,ind2);
if  isempty(V) %>
    Max = max(abs(diag(JD)))-1;
    LOCI= sort((diag(JD)))-1;
    sign_V =-1;
    vector = zeros(length(C),1);
else
    
    if V(C>0)
        V = V/V(C>0);
    end
    vector =equi_type*V;
    Max = MaxEig;
    LOCI= sort((diag(JD)))-1;
    sign_V =sign(C*A*vector);
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
if (m ~= n) || (n < 2)
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