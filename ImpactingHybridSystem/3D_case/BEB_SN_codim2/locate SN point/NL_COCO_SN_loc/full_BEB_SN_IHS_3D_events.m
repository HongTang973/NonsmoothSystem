function y = full_BEB_SN_IHS_3D_events(x, par, event)
%IMPACT_EVENTS   'hspo'-compatible encoding of event function.

switch event
  case 'impact'
    y = x(1,:);
  case 'zero_v'
    y = [1,0,0]*full_BEB_SN_IHS_3D_ode(x, par, '');
end

end
