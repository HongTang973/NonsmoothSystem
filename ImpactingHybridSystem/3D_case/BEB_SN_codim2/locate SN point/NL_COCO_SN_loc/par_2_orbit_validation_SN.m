% check if this will really give the same LCO
clc
close all
%% define the system
global const
scale_1     = 1;%-0.1; 
scale_2     = 0.3;% 0.01; 
scale_3     = 0.3;% 0.01; 
% from the formulae, we conclude that the c22 c23 c33 should be zero to get
% uniformly linear impact restitution law

A00 = scale_1* [0.25    0.1     0.2;
                0       -0.2    0.15;
                0.15    0.2     -0.1];

% B00*[x1*x2; x2*x3; x3*x1]
B00 = scale_2* [-0.6438    0   -0.5821;
                -0.2807    -0.3283  0.8103;
                -0.8866   -0.6487   0.3508];
% C00*[x1^2; x2^2; x3^2]
C00 = scale_3* [0.25*4     0       0;
                -0.5    -0.3   0.35;
                0.45   -0.6   0.1];
const.A00   = A00;
const.B00   = B00;
const.C00   = C00;
%% 
par = [-0.1,0.2,-0.5,1.78192697901049,1.6,4.79596124049338,0.02,0.0955685096979647]';


refine = 1;
options = odeset('Events',@(t,y) full_BEB_SN_IHS_3D_events_ode45_compat(t,y),...
    'RelTol',1e-12,'AbsTol',1e-12,'Refine',refine);
f_odes  =@(t,y) full_BEB_SN_IHS_3D_ode(y, par, '');
%
t0 = 0;
tmax   = 200;
y0     = [0; 0.0770926271121134; 0.0607559985948500];
% start calculation
tout    = t0;

yout    = y0.';
teout   = [];
yeout   = [];
yeout0  = [];
ieout   = [];

while t0 < tmax 
    
    [t,y,te,ye,ie] = ode45(@(t,y) f_odes(t,y),[t0 tmax],y0,options);
     nt    = length(t);
  
    % output history
    tout  = [tout; t(2:nt)];
    yout  = [yout; real(y(2:nt,:))];
    
    % output events
    teout = [teout; te];
    yeout = [yeout; ye];
    ieout = [ieout; ie];
    if size(y0,1)>1
        yeout0 = [yeout0; y0'];
    else
        yeout0 = [yeout0; y0];
    end
    t0=t(end);
    x=y(end,:)';
    if abs(t0-tmax) < 1e3*eps
        break
    end
    y0 = full_BEB_SN_IHS_3D_resets(x, par, 'bounce');
    tout = [tout;t0 + 1/1e4];
    yout = [yout; y0'];
    
end

figure
subplot(121)
plot(tout, yout(:,1))
subplot(122)
plot(yout(:,1), yout(:,2))
