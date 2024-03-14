function dRdP = Shilnikov_IHS_resets_dp(x, par, reset)
%IMPACT_RESETS   'hspo'-compatible encoding of reset function.
[~,B,C,~]           =   par2NForm_Shilnikov_DummyVar(par);
DBDP                =  [0, 0, 0, 0, 0, 0, 0, 0;
                        0, 0, 0, 0, 0, 0, 1, 0;
                        0, 0, 0, 0, 0, 0, 0, -1];
%
dRdP                =  - ( B*C*Shilnikov_IHS_DFDP(x, par, '') + ...
    C*Shilnikov_IHS_ode(x, par, '')*DBDP);
end
