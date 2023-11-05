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
function IC= IC_generator(T,R,A,C,equi_type)
method  = true;
%> monodromy matrix construction
Mono_A = expm(A*T);
%> the determinant condition
M = R*Mono_A- eye(length(C));
%> the straight way to calculate the determinant of the characteristic
%> matrix
f_1 = det(M);

PA = R*Mono_A;


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
% [ind, v]=find(abs(eig_A-1)<1e-9);
[v, ind]=min(abs(eig_A-1));


%
if ~isempty(ind)
    eig_A = eig_A(ind);
    JV = real(JV(:,ind));
    %>
else
    JV = real(JV(:,abs(eig_A-1)<1e-3));
end
%>
if size(JV,2) ~=1
    if norm(JV(:,1)-JV(:,2))<1e-6
        JV = JV(:,1);
    else
        error('There are two vectors corresponding to the unit multiplier!')
    end
end
% selelct real eigenvector
% JV = JV(:,sum(JV==real(JV))==length(C));
if JV(C>0)
    IC = JV./JV(C>0);
    IC = real(IC*equi_type);
end