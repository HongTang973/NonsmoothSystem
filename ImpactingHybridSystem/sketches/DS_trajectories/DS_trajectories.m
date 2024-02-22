%> this script use the 3D case code to show the trajectory interacting with
%> the discountinuity surface several times
close all
equi_type         = 1;
point_1           = [-0.1,0.2,-0.525,1.78192697901049+0.0955685096979647,1.6,4.7959612404933]';
[A,B,C,R, T_2_det]  = par2NForm_Lienard(point_1);
x_00                = IC_generator(point_1(6),R,A,C,equi_type)
prob                = BEB_SN_IHS_3D_par2prob(point_1,equi_type);
%
prob.IC             = 1.2*x_00;
prob.tspan          = [1 30];
[tout,yout,yeout0,teout,yeout,ieout]= Single_DS_IHS_INTEGRATION(prob);

% figure
% plot(tout, yout(:,3))

%> 

te = teout(ieout==1);
FIG1 = figure; hold on
ind0  = 1;
for i = 1:length(te)
    te_ind = tout>=te(i) & tout<=te(i);
    ie = find(te_ind>0);
    plot3(yout(ind0:ie,1),yout(ind0:ie,2),yout(ind0:ie,3),'k-','LineWidth',1.4)
    ind0 = ie + 1;
end
plot3(yout(ind0:end,1),yout(ind0:end,2),yout(ind0:end,3),'k-','LineWidth',1.4)

box off
set(gca,'XTick',[],'YTick',[],'ZTick',[])
xlim([0.5 2])
% patch
box0_x = [1 1 1 1 1];
box0_y = [-0.5 2.5 2.5 -0.5 -0.5 ];
box0_z = [-0.1 -0.1 0.3 0.3 -0.1 ];

V_0 = [box0_x(1), box0_y(1),box0_z(1);
    box0_x(2), box0_y(2),box0_z(2);
    box0_x(3), box0_y(3),box0_z(3);
    box0_x(4), box0_y(4),box0_z(4)];
F_0 = [1 2 3 4];

h2_3 =patch( 'Faces',F_0, 'Vertices',V_0,'EdgeColor',[1 1 1],'Linewidth',2, 'FaceColor', [0.1,0.1,0.3], 'FaceAlpha', 0.2,'displayname','$\Sigma$');
view(64.2,-48)
