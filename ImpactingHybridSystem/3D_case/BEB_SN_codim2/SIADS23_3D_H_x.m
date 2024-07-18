%>
function [value,isterminal,direction] = SIADS23_3D_H_x(t,y,prob)
value       =  real(y(1) - prob.equi_type);                %  detect switching point
isterminal  =  1;                              %  stop the integration
direction   = -1;
end