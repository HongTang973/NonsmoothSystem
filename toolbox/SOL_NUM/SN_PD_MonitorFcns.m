%> monitor the multipliers for the period one cycles in impacting hybrid
%> systems to check the stability and the PD/SN bifurcation
function [Event_fvals, S_US, p1_det ]= SN_PD_MonitorFcns(prob, par)
% added the return value det in the third place --> hopefully won't
% conflict with the previous code, since the previous codes don't require
% the third return var  ---- 3/9/224 Peter
[A,B,C,R]           = prob.par2prob(par);
T                   = par(6);
x_00                = IC_generator(T,R,A,C,1);
%> compute the Floquet Multipliers

F_ls                = @(y) A*y;
x_0                 = inv(R)*x_00;
%> consider the correction of the saltation matrix
correction  = ((F_ls(x_00)-R*F_ls(x_0))*C)/(C*F_ls(x_0));
P           =   R + correction;
p1_det      =   det(R*expm(A*T) - eye(length(C)));
PA          =   P*expm(A*T);
if any(any(isnan(PA))) || any(any(~isfinite(PA)))
    keyboard
end
[V2,D2]     =   eig(PA);
Salt_p      =   diag(D2);
%>
f_p         =   A*x_00;
ind         =   abs(Salt_p-1)<1e-6; %> get the trival unit multiplier
tmp         =   1:1:length(C);
index       =   tmp(ind);
%> locate the trivial unit multiplier
P_ind         = [];
cos_ang       = real(f_p'*V2 ./(vecnorm(V2)*norm(f_p)));
[~,P_ind]     = min(abs(abs(cos_ang) -1));
% if ~isempty(index) %> the perfect case when the solution is with det =0
%     for  i = 1:length(index)
%         % modified 2/12/2023 the cos(angle) should be real
%         v_i = V2(:,index(i));
%         cos_ang = real(f_p'*v_i/(norm(v_i)*norm(f_p)));
%         if abs(cos_ang^2-1)<1e-9 %> this means the multiplier is the trival one
%             P_ind            = index(i);
%             %< modified by Peter: SN point just need one trivial multiplier
%             break;
%         end
%     end
% else %> for the case when the perturbation on the parameter around the det =0 giving det= SN
%     cos_ang       = real(f_p'*V2 ./(vecnorm(V2)*norm(f_p)));
%     [~,P_ind]     = min(abs(abs(cos_ang) -1));
% end
%>
if isempty(P_ind)
    keyboard
end
Salt_p(P_ind) = 100; %> set the trival multiplier to a dummy value for
%> ease of the continuation
ind_trivial_multiplier = P_ind;
%> added by peter 4/11/2023 to check that if the cycle is stable and return the flag
S_US_ctr_err   = 1e-6;
S_US_indicator = Salt_p;
S_US_indicator(ind_trivial_multiplier) = [];
%> debug
% if any(abs(S_US_indicator)> 10)
%     keyboard
% end
%> determine the stability of the LCO
if any(abs(S_US_indicator)> 1 + S_US_ctr_err )
    S_US = 'US';
else
    S_US = 'S';
end
%> in case for some imaginary multipliers
% for j = 1:length(Salt_p)
%     if ~isreal(Salt_p(j))
%         Salt_p(j) = 100*sign(Salt_p(j));
%     end
% end
Salt_p              = sort(Salt_p);
%> modified by peter 3/11/2023 -- sove the problem that there is some cases
%> where the complex multipliers' real part cross the -1 or 1
% Event_fvals         = real([Salt_p-1; Salt_p+1]) + abs(imag([Salt_p-1; Salt_p+1]));
%> 30/11/2023 modified to resort
Event_fvals         = [sort(real(Salt_p-1) + abs(imag(Salt_p-1)));...
    sort( real(Salt_p+1) + abs(imag(Salt_p+1)) ) ];

% with this algorithm, the fake crossings will never be detected when the
% imaginary part is nonzero
% Event_fvals         =[Salt_p-1; Salt_p+1];

%> event function added on 05/12/2023 to detect the disappearance of the LCO
% with zero velocity
% Event_fvals  = [Event_fvals; C*A*x_00];
end