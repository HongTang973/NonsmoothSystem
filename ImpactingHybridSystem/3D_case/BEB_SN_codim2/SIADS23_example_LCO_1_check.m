%% > for given first point of LCO
clc
close all
equi_type   = 1;
fs          = 40;
%% stability check

% T    = roots(1);

[A,B,C,R]  = par2NForm(par_deli);
T          = par_deli(6);
F_ls = @(y) A*y;
% T = 5.847542373676547;
fprintf('For the first LCO: \n')
[sign_V,LOCI,Max,x_00,F_1,COND] = LCO_Det_search(T,R,A,C,equi_type);
% x_00 = 2*x_00;
x_0                 = inv(R)*x_00;
% C*A*x_0
% C*A*x_00
correction          = ((F_ls(x_00)-R*F_ls(x_0))*C)/(C*F_ls(x_0));
P                   = R + correction;
Mono_A              = expm(A*T);

%> 
% P2                  = R*(eye(3) - F_ls(x_00)*C/(C*F_ls(x_00)));
% PA2                 = P2*Mono_A;


%> 
PA      = P*Mono_A;
RA      = R*Mono_A;

[V1,D1] = eig(RA);
Mono_p  = diag(D1)
V1

[V2,D2] = eig(PA);
[V2_L,D2_L] = eig(PA');
Salt_p  = diag(D2)
Salt_p_L  = diag(D2_L)
V2
% V2_L
v = real( V2(:,2) )
w = real( V2_L(:,2) )
fprintf('The dot product between the left and right eigenvector is %g. \n', w'*v)
% V2_L'*V2
% w'*PA*v

%> 

fp = A*x_00;
fprintf('Check the A*y is the eigenvector of the Phi_x of unit eigenvalue \n')
v'*fp/norm(fp)/norm(v)
% w = w/(v'*w)

% eig(PA2)

% subV2 = V2(:,1:2)';
% 
% dim2_jac = subV2*PA*pinv(subV2)
% [V_dim2, D_dim2] = eig(dim2_jac)
% 
% B2 = V_dim2'*subV2

%> 


% eigV_fold = V2(:,2)';
% fold_jac = eigV_fold*PA*pinv(eigV_fold)


% sub_RA = RA(2:3, 2:3);
% eig(sub_RA)
% 
%> project the map on to the poincare section
Kel = vec2kernel(real(V2(:,2)))
TT = Kel'*PA*Kel
[V_T, D_T]= eig(TT)
V_T(:,1)'*V_T(:,2)
% sub_PA = PA(2:3, 2:3);
% eig(sub_PA)

%> numerical check the return 
F_ls_ode = @(t,y) A*y;
InitCond = x_00;

% C*A*InitCond

refine = 1;
options = odeset('Events',@(t,y) H_x(t,y),...
        'RelTol',1e-12,'AbsTol',1e-12,'Refine',refine);
% [tout,yout,yeout0,teout,yeout,ieout]=...
%      Single_DS_Impacting_Hybrid_system_integration(A,R,C,InitCond,tspan,fs,equi_type);
%  figure
% plot(tout, yout(:,1))

[t1,y1,te,ye,~] = ode45(@(t,y) F_ls_ode(t,y), [0 10*T], InitCond, options);
%>
figure
plot(t1, y1(:,1))

%% > use the iteration to numerically evaluate the Poincare map


% -------------------------  the second LCO ---------------------------- %
%>
% T = roots(2);
% fprintf('For the second LCO: \n')
% [sign_V,LOCI,Max,x_00,F_1,COND] = LCO_Det_search(T,R,A,C, equi_type);
% x_0 = inv(R)*x_00;
% correction = ((F_ls(x_00)-R*F_ls(x_0))*C)/(C*F_ls(x_0));
% P =R + correction;
% Mono_A = expm(A*T);
% PA      = P*Mono_A;
% [V1,D1]=eig(R*Mono_A);
% Mono_p =diag(D1)
% % abs(Mono_p)
% [V2,D2]=eig(PA);
% Salt_p =diag(D2)

% --- calculate the Poincare map ---
 RMap = @(y0) R*y0;
 
[dTdy2, dTdy3] = numerical_dTdy(x_00, T, fs, F_ls_ode, RMap, @H_x)
 % ------  linear part of the 2D ------- 
 [J_pMP,eigD_pMP, eigV_pMP] = Numerical_2D_PMP(x_00, T, fs, F_ls_ode, RMap, @H_x)
 w = eigV_pMP(:,2)
 eigV_pMP(:,2)'*eigV_pMP(:,1)
%  w = eigV_pMP(:,1)
 
 [lambda, c] = Numerical_2D_PMP_normal_form(w,x_00, T, fs, F_ls_ode, RMap,  @H_x)

 [J_pMP,eigD_pMP, eigV_pMP] = Numerical_3D_PMP(x_00, T, fs, F_ls_ode, RMap, @H_x_mu)
 PA
 RA
 
 % ------  nonlinear part  ------- 
 
 
 function [value,isterminal,direction] = H_x(t,y)
 value = real(y(1) - 1);                %  detect switching point
 isterminal = 1;                              %  stop the integration
 direction  = -1;
 end
 
function [value,isterminal,direction] = H_x_0(t,y)
 value = real(y(1) - 1);                %  detect switching point
 isterminal = 0;                              %  stop the integration
 direction  = -1;
end

%>
function [value,isterminal,direction] = H_x_mu(t,y)
 value = real(y(1) - mu);                %  detect switching point
 isterminal = 0;                              %  stop the integration
 direction  = -1;
end