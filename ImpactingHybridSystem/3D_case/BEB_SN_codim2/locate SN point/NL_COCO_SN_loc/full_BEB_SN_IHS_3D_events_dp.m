function y = Shilnikov_IHS_events_dp(x, par, event)
%IMPACT_EVENTS   'hspo'-compatible encoding of event function.

switch event
  case 'impact'
    y = zeros(1,8);
  case 'zero_v'
    y = [1,0,0]*Shilnikov_IHS_DFDP(x, par, '');
end

end
