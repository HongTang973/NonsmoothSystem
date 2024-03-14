function dRdX = Shilnikov_IHS_resets_dx(x, par, reset)
%IMPACT_RESETS   'hspo'-compatible encoding of reset function.
[~,B,C,~]       =   par2NForm_Shilnikov_DummyVar(par);
    dRdX           =   eye(3) - B*C*Shilnikov_IHS_DFDX(x, par, ''); 
end
