function [data, res] = full_BEB_SN_IHS_3D_add_IP(prob, data, command, varargin)
%DUFF_ADD_IP   Slot function: Add to bifurcation data.
%
% Store initial point on first trajectory segment to bifurcation data

res = {};
switch command
  case 'init'
    res   = 'X0';
  case 'data'
    chart = varargin{1};
    [data, uidx] = coco_get_func_data(prob, 'hspo.orb.bvp.seg1.coll', ...
      'data', 'uidx');
    res  = chart.x(uidx(data.coll_seg.maps.x0_idx));
end

end
