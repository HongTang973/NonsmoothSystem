function f = BilinearFunction(B,x,y)
dim_sys = size(B,1);
if dim_sys~= length(x) || dim_sys~= length(y); error('Dismatch in size!');end
%> 
% krov = [];
% for i = 1:dim_sys; krov = [krov; x(i).*y]; end
krov = zeros(dim_sys^2,1);
for i = 1:dim_sys
    for j = 1:dim_sys
        krov((i-1)*dim_sys + j) = x(i)*y(j);
    end
end
%
f = B*krov;
end