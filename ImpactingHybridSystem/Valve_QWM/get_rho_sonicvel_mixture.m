function [rho_m, a_m, alpha_g] =...
    get_rho_sonicvel_mixture(xg,rho_g,rho_l, a_g, a_l)


alpha_g = xg./(xg + (1 - xg)*rho_g/rho_l);

rho_m   = alpha_g*rho_g + (1- alpha_g)*rho_l;

temp_1_a2 = rho_m.*(alpha_g/rho_g/a_g^2 + (1-alpha_g)/rho_l/a_l^2);

a_m = sqrt(1./temp_1_a2);
end