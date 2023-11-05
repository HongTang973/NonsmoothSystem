function F = SIADS23_SN_zero_Fcns( par)

[A,B,C,R, T_2_det,T_2_d_det_dt]       = par2NForm_Lienard(par);
F               = [T_2_det(par(6) ); T_2_d_det_dt(par(6) )];

end