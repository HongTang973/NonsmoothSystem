function y = Shilnikov_IHS_events_dx(x, par, event)
%IMPACT_EVENTS   'hspo'-compatible encoding of event function.

switch event
  case 'impact'
    y = [1 0 0];
  case 'zero_v'
    y = [1,0,0]*Shilnikov_IHS_DFDX(x, par, '');
end

end
