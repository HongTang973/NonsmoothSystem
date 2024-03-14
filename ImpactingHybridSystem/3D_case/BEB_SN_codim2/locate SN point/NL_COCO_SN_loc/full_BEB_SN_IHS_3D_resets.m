function y = full_BEB_SN_IHS_3D_resets(x, par, reset)
%IMPACT_RESETS   'hspo'-compatible encoding of reset function.

[A,B,C,R]       = par2NForm_DummyVar(par);
switch reset
    case 'bounce'
        y       = x - B*C*full_BEB_SN_IHS_3D_ode(x, par, '');
    case 'zero_v'
        y       = x;
        % no action
end

end
