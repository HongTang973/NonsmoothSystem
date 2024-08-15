function [Event_fvals, S_US ]= SIAMDS24_SN_codim_MonitorFcns(prob, par)
        [Event_fvals, S_US ]= SN_PD_MonitorFcns(prob, par);
%         for further codim2 analysis : the line zeta = 0.64 is chosen
        Event_fvals        = [Event_fvals(4:6); par(7)-0.025;par(7)-0.02; par(4)-1.90];
        % Event_fvals        = [Event_fvals(1:3); par(4)-0.653363855944258-0.02;100; 100];

end