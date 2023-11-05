%>
function [value,isterminal,direction] = SIADS23_3D_H_x(t,y)
value       =  real(y(1) - 1);                %  detect switching point
isterminal  =  1;                              %  stop the integration
direction   = -1;
end