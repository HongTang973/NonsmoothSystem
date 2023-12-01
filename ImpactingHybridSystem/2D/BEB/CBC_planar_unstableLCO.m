% CBC technique to show the unstable LCO in planar impacting system
clc
clear
close all
% define the linear part
% As shown in the Exapmle 2.1, 
% xi = 1.05; r = 1.5;
% A  =  [0 1;-1 -2*xi];
a =0.5;b=-1;
A = [0,1;b,a];
r= 0.5;
MW = [0 0;0 -(1+r)];
R  = eye(size(A,1)) + MW;
C  = [1,0];

% define the useful function
[V1,D1]=eig(A)
exp(-real(D1(1,1))*pi/imag(D1(1,1)))
r*exp(real(D1(1,1))*pi/imag(D1(1,1)))
EA = @(T) real(V1*diag(exp(diag(D1)*T))*inv(V1));
sign_V = @(Y) sign(C*A*Y);
% determine the sampling frequency 
Omega = max(abs(diag(D1)));
fs = 100*ceil(2*(Omega/2/pi));
t0 =100;
tspan = [0 t0];
% DNS   
% if equi_type == 'admissible'
%         pre = -1;
%  elseif equi_type == 'pseudo'
%         pre = 1;
%  end
equi_type =   1;                                                                           
%
InitCond = [equi_type 10]';
discrepancy = 1;
iter =0;
while discrepancy > 1e-4 && iter <20
[tout,yout,yeout0,teout,yeout,ieout]=...
    Single_DS_Impacting_Hybrid_system_integration_CBC(A,R,C,InitCond,tspan,fs,equi_type);
discrepancy =norm(yeout0(end,:)' -InitCond)/norm(InitCond);
% y0_ = InitCond;
% if (length(teout)>1)
%     y0_ = y0_+R*(EA(teout(end)-teout(end-1))*y0_-yeout(end,:)');
%     else
%     y0_ = y0_+R*(EA(teout(end)-0)*y0_-yeout(end,:)');    
% end
%     y0_(C>0) =equi_type;
InitCond=yeout0(end,:)';
iter = iter +1
end
%

figure
plot(tout,yout(:,C>0),'r-','linewidth',1.4)
title(['LCO ','''s stability'])
set(gca,'fontname','times new roman','fontsize',12)
xlabel('t/s')
figure
plot(yout(:,1),yout(:,2),'k-','linewidth',1.2)
title(['LCO ','''s phase portrait'])
set(gca,'fontname','times new roman','fontsize',12)
grid on