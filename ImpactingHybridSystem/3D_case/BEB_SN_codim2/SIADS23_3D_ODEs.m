%> whenever you define the function to form the flow operator and the
%> imapct map, you already define the index of the parameters, the
%> dimension of the system

function  Flow_ode =SIADS23_3D_ODEs(par)
%>
[A,B,C,R,T_2_det] = par2NForm_Lienard(par);
Flow_ode         = @(t,y) A*y;
%>
end