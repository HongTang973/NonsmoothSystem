% return the Floque multipliers of the monodary one and the one with
% saltation matrix
% 
function [Mono_p,Salt_p]=IC2Floque_Multipliers(T,x_00,R,A,C)

F_ls = @(y) A*y;
x_0 = inv(R)*x_00;
%> consider the correction of the saltation matrix
correction = ((F_ls(x_00)-R*F_ls(x_0))*C)/(C*F_ls(x_0));
P =R + correction;
Mono_A = R*expm(A*T);
PA = P*expm(A*T);
[V1,D1]=eig(Mono_A);
Mono_p =diag(D1);
[V2,D2]=eig(PA);
Salt_p =diag(D2);
end