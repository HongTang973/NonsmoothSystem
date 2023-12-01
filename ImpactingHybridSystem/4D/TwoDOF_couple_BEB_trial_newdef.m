% 2Dof coupled impact oscillator sweep
addpath ../../public_scripts
%*************** start parallel computation ***************%
% find the numbers of this 
% Initialization_for_Par;
% case 1: sweep the coupling parameters for non-bif case
% case 2: sweep the coupling parameters for bif case
clc
clear 
close all
% Define coupling effect parameter list
up_limit  = 1;
bot_limit = 0;
delta_v   = 0.001;
epsilon_list=bot_limit:delta_v :up_limit;

% Parameter Definition

%
T=300;
tspan=[0 T];fs=1000;
zeta1=-0.1; zeta2=0.1;
mu=1/3;   % w2/w1
eta=1/10;   % m2/m1
epsilon=0.5; % the perturbation
res=0.3;
omega=1;

m12=1;m21=m12*eta; m_c=epsilon*m21;    
c12=-2*zeta1*mu;c21=c12*eta;
k21=-mu^2;k12=k21*eta;
X1_0=0;dX1_0=0.5062;X2_0=0;dX2_0=0;
InitCond=[X1_0;X2_0;dX1_0;dX2_0];
% get the eigenvalue varying process

% time series
a1=-mu*zeta1;a2=-zeta2;b1=mu*sqrt(1-zeta1^2);b2=sqrt(1-zeta2^2);
A1=-(m12*(a2^2-b2^2)+c12*a2+k12);A2=2*a2*b2*m12+b2*c12;
B1=-(m21*(a1^2-b1^2)+c21*a1+k21);B2=2*a1*b1*m21+b1*c21;
phi_1=atan((a1*X1_0-dX1_0)/b1/X1_0);
phi_2=atan((a2*X2_0-dX2_0)/b2/X2_0);%


% derived compound parameters
if X1_0~=0
    Am10=X1_0/cos(phi_1);
else
    Am10=dX1_0/(a1*cos(phi_1)-b1*sin(phi_1));
end
Am11=0.5*Amplitude(a2-a1,b2-b1,A1,-A2)/b1;
Am12=0.5*Amplitude(a2-a1,b2-b1,A2,A1)/b1;
Am13=0.5*Amplitude(a2-a1,b2+b1,A1,-A2)/b1;
Am14=0.5*Amplitude(a2-a1,b2+b1,A2,A1)/b1;

if X2_0~=0
    Am20=X2_0/cos(phi_2);
else
    Am20=dX2_0/(a2*cos(phi_2)-b2*sin(phi_2));
end
Am21=0.5*Amplitude(a2-a1,b2-b1,-B1,B2)/b2;
Am22=0.5*Amplitude(a2-a1,b2-b1,B2,B1)/b2;
Am23=0.5*Amplitude(a2-a1,b2+b1,-B1,-B2)/b2;
Am24=0.5*Amplitude(a2-a1,b2+b1,-B2,B1)/b2;

% get the matrix 
% sweep the parameters

MHH=[1,0;0,1]+epsilon*[0 m12;m21,0];
K=[mu^2,0;0,1]+epsilon*[0 k12;k21,0];
C=[mu*2*zeta1,0;0,2*zeta2]+epsilon*[0 c12;c21,0];
% eig(K,M)
%
Preload=[0;0;omega;0];
%
A_matrix=[zeros(2,2),eye(2);-MHH\K,-MHH\C];
A_matrix\Preload
eig(A_matrix)
% sweep the list 

% time integration
[flag,tout,yout,yeout0,teout,yeout,ieout]=...
    DirectInt(m_c,res,A_matrix,InitCond,tspan,MHH,Preload);
%Period
Period=teout(2:end)-teout(1:end-1);
%
t_kept=T-100;
abstract_ind=(teout>t_kept);
if ~isempty(abstract_ind) 

end
%
figure
subplot(211)
plot(tout,yout(:,1),'k-','linewidth',1.2)
% hold on
% plot(tout,X20(t),'r','linewidth',1.2)
legend('Mass1-num')
subplot(212)
plot(tout,yout(:,2),'k-','linewidth',1.2)
% hold on
% plot(tout,X20(t),'r','linewidth',1.2)
legend('Mass2-num')

figure
subplot(211)
plot(tout,yout(:,3),'k-','linewidth',1.2)
% hold on
% plot(tout,X20(t),'r','linewidth',1.2)
legend('Mass1-num')
subplot(212)
plot(tout,yout(:,4),'k-','linewidth',1.2)
% hold on
% plot(tout,X20(t),'r','linewidth',1.2)
legend('Mass2-num')
% 
figure(7)
subplot(211)
plot(yout(end-50*fs:end,1),yout(end-50*fs:end,3),'k-','linewidth',1)
 hold on
% plot(tout,X20(t),'r','linewidth',1.2)
legend('Mass1-num')
legend('Mass1-num')
subplot(212)
plot(yout(end-50*fs:end,2),yout(end-50*fs:end,4),'k.-','linewidth',1)
hold on
% plot(tout,X20(t),'r','linewidth',1.2)
legend('Mass2-num')



% semi analytic results and get contarst with DN results
V_LCO_1=1.515804072208207;
T01=14.257720657037709;
X0=0;
dX0=V_LCO_1;
Y0=-0.010003410173754;
dY0=-0.066272212711759;
InitCond=omega*[X0;Y0;dX0;dY0];
vect_preload=[0;0;omega;0];
[V,D]=eig(A_matrix);
mcoef=V*diag(V\(InitCond+A_matrix\vect_preload));
eigen_value=diag(D);
%ini=A_matrix\vect_preload(:,1)
% define beta & dot_beta time history function
beta     =@(mcoef,eigen_value,t) mcoef(1,:)*exp(eigen_value*t)-sum(mcoef(1,:))+X0;
dot_beta =@(mcoef,eigen_value,t) mcoef(1,:)*(eigen_value.*exp(eigen_value*t));
Y    =@(mcoef,eigen_value,t) mcoef(2,:)*exp(eigen_value*t)-sum(mcoef(2,:))+omega*Y0;
dY =@(mcoef,eigen_value,t) mcoef(2,:)*(eigen_value.*exp(eigen_value*t));
T_list=0:1/fs:T01;
DisX=beta(mcoef,eigen_value,T_list);
VelX=dot_beta(mcoef,eigen_value,T_list);
DisY=Y(mcoef,eigen_value,T_list);
VelY=dY(mcoef,eigen_value,T_list);
figure(7)
subplot(211)
plot(DisX,VelX,'r-','displayname','Analytic');
subplot(212)
plot(DisY,VelY,'r-','displayname','Analytic')
figure
plot(T_list,DisX)


n=10;
figure
subplot(211)
plot(yout(1:n*fs,1),yout(1:n*fs,3),'k-','linewidth',1)
 hold on
% plot(tout,X20(t),'r','linewidth',1.2)
% legend('Mass1-num')
legend('Mass1-num')
subplot(212)
plot(yout(1:n*fs,2),yout(1:n*fs,4),'k.-','linewidth',1)
hold on
% plot(tout,X20(t),'r','linewidth',1.2)
legend('Mass2-num')










% sub_functions
function [flag,tout,yout,yeout0,teout1,yeout1,ieout]=DirectInt(m_c,res,A_matrix,InitCond,tspan,MHH,Preload)
flag=0;
fs=1000;
%
tstart=tspan(1);
tfinal=tspan(2);
%
y0=InitCond;

% check if start at the boundary
if abs(y0(1))<1e-3 && y0(3)<0
    y0=impact_map(MHH,res,y0);
end
%
opt_ls=odeset('Events',@(t,y) Detect_impact(t,y),'RelTol',1e-12,'AbsTol',1e-12,'InitialStep',1e-4,'maxstep',1e-3);
opt_stuck=odeset('Events',@(t,y) Detect_unstickevent(t,y,A_matrix,Preload),'RelTol',1e-6,'AbsTol',1e-6,'InitialStep',1e-4);
%

% function during region I
f_ls= @(t,y,A_matrix,Preload) A_matrix*y+Preload;
hx_= @(t,y) [0,0,1,0];  % H_xx=0; we gwt a(x)=HxFx*F
Wx_= @(t,y) (1+res)*[0 0 -1 m_c]';
f_sk= @(t,y,A_matrix,Preload) f_ls(t,y,A_matrix,Preload)+ (hx_(t,y)*f_ls(t,y,A_matrix,Preload))*Wx_(t,y)/(1+res)  ;

% time integration
tout = tstart;
yout = y0.';
yeout0 = []; % record the initial v for the impact
teout1 = [];
yeout1 = [];
ieout =  [];

BeStuck=false;

while tout(end)<tfinal
    % define out put time serias
    try
        timespan=[tstart:1/fs:tfinal];
    catch ME
        keyboard
    end
    %     timespan=[tstart tfinal];
    
    % exit case I
    if length(timespan)<2
        % last time step
        break;
    end
    if BeStuck==false
    % Solve until the first terminal event: hit the boundary.
    [t,y,te,ye,ie] = ode45(@(t,y) f_ls(t,y,A_matrix,Preload),...
        timespan,y0,opt_ls);
    
    else % if it gets stuck motion at this boundary 
        f_ls(t,y0,A_matrix,Preload)
        f_sk(t,y0,A_matrix,Preload)
     [t,y,te,ye,ie] = ode45(@(t,y) f_sk(t,y,A_matrix,Preload),...
        timespan,y0,opt_stuck);
   
    end
    
    % debug
%     figure(1)
%     plot(t,y(:,3))
%     figure(2)
%     plot(t,y(:,1))
    % Accumulate output.  This could be passed out as output arguments.
    
    nt    = length(t);
    %     if t(end)>2
    %         keyboard
    %     end
    
    % output history
    tout  = [tout; t(2:nt)];
    yout  = [yout; y(2:nt,:)];
    % stop when divergence
    if any(abs(yout(end,:))>100)
        flag=1;
        break
    end
    
    % output events
    teout1 = [teout1; te];
    yeout1 = [yeout1; ye];
    if size(y0,1)>1
    yeout0 = [yeout0; y0'];
    else
        yeout0 = [yeout0; y0];
    end
    ieout = [ieout; ie];
    %evaluate new initial condition for next integration step if possible
    if isempty(ie)
        % means that no event detected and the integration finished till
        % the time destination, so just end the intgration and return
        % result to the mainfunction
        break
    elseif BeStuck==true
        % if previous integration is sliding, then it shows the sliding
        % motion stopped and continue intefration in f_ls
        y0=y(nt,:);
        tstart=t(nt);
        BeStuck=false; % stuck state over
        continue % return to the biginning of this time integration loop
    end
    
   %  ------ else applay impact map ---------%
   t0=t(end);
   x=y(end,:)';
   vx=x(3);
   Fx=f_ls(t0,x,A_matrix,Preload);
   ax=[0 0 1 0]*Fx;
   if vx>0 && t0<tfinal
       % this case is abnormal and need further attention 
   else % negative normal impact velocity
     
   end
   % end this chattering sequence when the velocity meets the threshold
   
   if abs(vx)<1e-2 && ax<0 &&t0< tfinal
       % apply chttering mapping 
       Wx=(1+res)*[0;0;-1;m_c];
       xstar=x+(1/(1-res) * ( 2 * Fx * res/ax + Wx ))*vx;
       tstar =t0+1/(1-res) * ( 2*res/ax ) * vx;
       if abs(xstar(1))<0
            % this map occurrs when hitting the boundary but it's wrong
            % when the flag freedom penetrates the boundary, so we need fix
            % this mannuly with physical meaning
            
            %             disp('wrong chattering map')
            xstar(1)=0;
            xstar(3)=0;
       end
        % after the chattering sequence ends
   BeStuck=true;
   % append the map resulted state to the output
   tout=[tout;tstar];
   yout=[yout;xstar'];
   y0=xstar;
   else
       % normal impact 
        y0=y(nt,:);
        y0=impact_map(MHH,res,y0);
        
        
        tstart=t(nt)+1/fs;
        yout=[yout;y0];
        tout=[tout;tstart];
        BeStuck =false;
   end
  
  %   
end

% to detect the contact
    function [value,isterminal,direction]=Detect_impact(t,y)
        value = y(1);       %  detect event
        isterminal = 1;        %  stop the integration
        direction  = -1;
    end
% to detect the unsitck event
    function [value,isterminal,direction]=Detect_unstickevent(t,y,A_matrix,Preload)
        F=A_matrix*y+Preload;
        value = [0 0 1 0]*F;      %  detect event when the acc becomes positive
        isterminal = 1;           %  stop the integration
        direction  = 1;
    end
%
    function y0=impact_map(MHH,res,y0)
        delta_v=MHH(2,1)/MHH(2,2)*(1+res)*y0(3);
        y0(3)=-res.*y0(3);
        y0(4)=y0(4)+delta_v(1);
    end
%
end


% subfunction 2
function Am=Amplitude(A,B,C,D)
%
Am=(C*A+D*B)/(A^2+B^2);

end

