function [Event_fvals, S_US ]= SIAMDS24_SN_codim_MonitorFcns_Dummy(prob, par)
        [Event_fvals, S_US ]= SN_PD_MonitorFcns(prob, par);
%         for further codim2 analysis : the line b2 = 1.9 is chosen
        Event_fvals        = [Event_fvals(4:6); par(7)-0.025;par(7)-0.02; 1.90-par(4)-par(8)];
        

end