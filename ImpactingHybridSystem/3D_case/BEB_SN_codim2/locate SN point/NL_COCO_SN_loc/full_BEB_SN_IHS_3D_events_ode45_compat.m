function [value,isterminal,direction] =... 
    full_BEB_SN_IHS_3D_events_ode45_compat(t,y)
%IMPACT_EVENTS   'hspo'-compatible encoding of event function.

% Locate the time when flap hit through boundary in both directions and
% stop integration when isterminal ==1 Events: switch boundary; zero heave
% velocity; zero pitch velocity; zero flap velocity
value       = y(1);
%  detect switching point
isterminal  = 1;                              %  stop the integration
direction   = -1;                             %  mark with event funct-
%  ion adding direction

end
