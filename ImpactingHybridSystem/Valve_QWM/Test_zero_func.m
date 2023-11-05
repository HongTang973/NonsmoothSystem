T = 0.502000817471173;
par = [343.14, 0.78, 0.05 , 0.5,  T]';

[prob_det] = ZeroFunctions_Valve_QWM (par)

par = [343.14, 0.78, 0.05 , 1.058,  100, 1]';
[prob_det] = ZeroFunctions_Valve_QWM_HB(par)

par = [343.14, 0.78, 0.05 , 4.003,  100, 0]';
[prob_det] = ZeroFunctions_Valve_QWM_HB(par)


T = 0.562225056853600;
par = [343.14, 0.78, 0.05 , 0.5,  T]';
[prob_det] = ZeroFunctions_Valve_QWM (par)

xg_1= 0:1e-10:1e-8;
xg_2 = 1e-8 : 1e-8:1e-6;
xg_3 = 1e-6 : 1e-6 : 1e-4;
xg_4= 1e-4:1e-4:1;
xg = [xg_1, xg_2,xg_3, xg_4];
rho_g   = 1.2;
rho_l   = 1000;
a_g     = p.a;
a_l     = 1300;
[rho_m, a_m, alpha_g] = ...
    get_rho_sonicvel_mixture(xg,rho_g,rho_l, a_g, a_l);
close all
figure
yyaxis left 
semilogx(xg, a_m)
ylabel('a_m')
yyaxis right 
semilogx(xg, alpha_g)
ylabel('alpha_g')
grid on
xlabel('xg')