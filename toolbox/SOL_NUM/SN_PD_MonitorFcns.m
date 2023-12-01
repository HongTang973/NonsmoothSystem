%> monitor the multipliers for the period one cycles in impacting hybrid
%> systems to check the stability and the PD/SN bifurcation
function [Event_fvals, S_US ]= SN_PD_MonitorFcns(prob, par)
[A,B,C,R]           = prob.par2prob(par);
T                   = par(6);
x_00                = IC_generator(T,R,A,C,1);
%> compute the Floquet Multipliers

F_ls                = @(y) A*y;
x_0                 = inv(R)*x_00;
%> consider the correction of the saltation matrix
correction  = ((F_ls(x_00)-R*F_ls(x_0))*C)/(C*F_ls(x_0));
P           =   R + correction;
PA          =   P*expm(A*T);
if any(any(isnan(PA))) || any(any(~isfinite(PA)))
    keyboard
end
[V2,D2]     =   eig(PA);
Salt_p      =   diag(D2);
%>
f_p         = A*x_00;
ind         = abs(Salt_p-1)<1e-6; %> get the trival unit multiplier
tmp         = 1:1:length(C);
index       = tmp(ind);
%> locate the trivial unit multiplier
for  i = 1:length(index)
    v_i = V2(:,index(i));
    cos_ang = f_p'*v_i/(norm(v_i)*norm(f_p));
    if abs(cos_ang^2-1)<1e-9 %> this means the multiplier is the trival one
        Salt_p(index(i)) = 100; %> set the trival multiplier to a dummy value for
        %> ease of the continuation
        ind_trivial_multiplier = index(i);
        %> added by peter 4/11/2023 to check that if the cycle is stable and return the flag
        S_US_ctr_err   = 1e-6;
        S_US_indicator = Salt_p;
        S_US_indicator(ind_trivial_multiplier) = [];
        if any(abs(S_US_indicator)> 1 + S_US_ctr_err )
            S_US = 'US';
        else
            S_US = 'S';
        end
        %> determine the stability of the LCO
    end
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

end