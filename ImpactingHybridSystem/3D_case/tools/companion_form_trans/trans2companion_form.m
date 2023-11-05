function [A_c,B_c,OZ] = trans2companion_form(A,B,OC)
n_d = size(A,1);
% OC = zeros(n_d,1);
% OC(1) = 1;
% OC = [1;1;1];
% OC : observing vector
temp = eye(n_d);
for i = 1:n_d-1
    temp = temp*(A');
    OC = [temp*OC(:,end),OC];
end

z = OC(:,end)'/OC;
OZ = z';
temp = eye(n_d);
for i = 1:n_d-1
    temp = temp*A;
    OZ = [temp*OZ(:,end),OZ];
end

A_c = OZ\A*OZ;
B_c = OZ\B;

end