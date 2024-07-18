function f = TrilinearFunction(C,x,y,z)
dim_sys = size(C,1);
if dim_sys~= length(x) || dim_sys~= length(y) || dim_sys~= length(z); error('Dismatch in size!');end
%> 
% krov = [];
% for i = 1:dim_sys; krov = [krov; x(i).*y]; end
% %> 
% Krov = [];
% for k = 1:dim_sys; Krov = [Krov; z(k).*krov]; end
krov = zeros(dim_sys^3,1);
for i = 1:dim_sys
    for j = 1:dim_sys
        for k = 1 : dim_sys
        krov((i-1)*dim_sys + j + dim_sys^2*(k-1)) = x(i)*y(j)*z(k);
        end
    end
end
%
f = C*krov;
end