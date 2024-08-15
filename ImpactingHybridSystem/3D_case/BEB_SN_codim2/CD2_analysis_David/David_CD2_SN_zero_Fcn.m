function P1_SN = David_CD2_SN_zero_Fcn(prob, par)
[Event_fvals, ~ ,p1_det]   = SN_PD_MonitorFcns(prob, par);
SN_det              = 1e-3*det(diag(Event_fvals(1:2)));
% we get rid of the third element which is 99
P1_SN               = [p1_det;SN_det];
end