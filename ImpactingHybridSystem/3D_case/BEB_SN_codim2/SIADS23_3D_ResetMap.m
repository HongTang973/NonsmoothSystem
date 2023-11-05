function  ResetMap = SIADS23_3D_ResetMap(par)
%>
[A,B,C,R] = par2NForm_Lienard(par);
ResetMap        = @(y0) R*y0;
end