% 2Dof coupled impact oscillator sweep
addpath ../../public_scripts
%*************** start parallel computation ***************%
% find the numbers of this
Initialization_for_Par;
% case 1: sweep the coupling parameters for non-bif case
% case 2: sweep the coupling parameters for bif case
clc
clear
close all
% Define coupling effect parameter list
up_limit  = 0.12;
bot_limit = 0.0;
delta   = 0.001;
epsilon_list=bot_limit:delta :up_limit;
%
num=1;
x_b='$\rm{Coupling \quad term}(\epsilon)$';
% x_b='$f$';
% define par loop parameter list
init_indx  =[];
par_list =[];
for i=1:num
    par_list=[par_list,epsilon_list];
end
% initcond
for j=1:num
    init_indx=[init_indx,j*3*ones(1,length(epsilon_list))];
end
% init_indx=[init_indx,-init_indx];
% Parameter Definition

%
T=200;
tspan=[0 T];fs=20;
zeta1=-0.1; zeta2=0.1;
mu=3;   % w2/w1
eta=10;   % m2/m1
epsilon=0.175; % the perturbation
res=0.3;


m21=1;m12=m21*eta; 
c21=2*zeta2*mu;c12=c21*eta;
k21=mu^2;k12=k21*eta;
% X1_0=0;dX1_0=3;X2_0=0;dX2_0=0;
% InitCond=[X1_0;X2_0;dX1_0;dX2_0];

% time series
a1=-zeta1;a2=-zeta2*mu;b1=sqrt(1-zeta1^2);b2=mu*sqrt(1-zeta2^2);

% get the eigenvalue varying process
% get the matrix
omega=1;
eigenvalue=[];
for i=1:length(epsilon_list)
    epsilon=epsilon_list(i);
    MHH=[1,0;0,1]+epsilon*[0 m12;m21,0];
    K=[1,0;0,mu^2]+epsilon*[0 k12;k21,0];
    C=[2*zeta1,0;0,2*zeta2*mu]+epsilon*[0 c12;c21,0];
    % eig(K,M)
    %
    Preload=[0;0;omega;0];
    %
    A_matrix=[zeros(2,2),eye(2);-MHH\K,-MHH\C];
    [V,D]=eig(A_matrix);
    selection=diag(D);
%     if i>1
%         discrepancy=abs(eigenvalue(:,end)-selection.');
%         [value,ind]=min(discrepancy);
%         selection=selection(ind);
%     end
    eigenvalue=[eigenvalue,selection];
end
figure
plot(epsilon_list,real(eigenvalue(1,:)),'r.',epsilon_list,real(eigenvalue(3,:)),'k.','markersize',5)
hold on
plot(epsilon_list,real(eigenvalue(2,:)),'r.',epsilon_list,real(eigenvalue(4,:)),'k.','markersize',5)
xlabel('$\rm{Coupling \quad term}(\epsilon)$','interpreter','latex')
ylabel('$\rm{Real \quad part}$','interpreter','latex')
plot(epsilon_list,0*real(eigenvalue(1,:)),'b.')
legend('$Rel(\lambda_1)$','$Rel(\lambda_2)$','interpreter','latex')
% ylim([-1 1])
set(gca,'fontname','times new roman')
figure
plot(epsilon_list,imag(eigenvalue(1,:)),'r.',epsilon_list,imag(eigenvalue(2,:)),'r.',epsilon_list,imag(eigenvalue(3,:)),'k.',epsilon_list,imag(eigenvalue(4,:)),'k.','markersize',5)% derived parameters



% get the matrix
% sweep the parameters
m1_v=[];
m2_v=[];
length_m1_v=[];
lenght_m2_v=[];
tic
parfor j=1:length(par_list)
%     try
    m1_v_sec=[];
    m2_v_sec=[];
    
    bif_para=par_list(j);
    epsilon=bif_para;
    omega=1;
    
%     bif_para=par_list(j);
%     omega=bif_para;
%     epsilon=0;

    m_c=epsilon*m21;
    
    disp(bif_para)
    X1_0=0;dX1_0=init_indx(j);X2_0=0;dX2_0=1;
    InitCond=[X1_0;X2_0;dX1_0;dX2_0];
    MHH=[1,0;0,1]+epsilon*[0 m12;m21,0];
    K=[1,0;0,mu^2]+epsilon*[0 k12;k21,0];
    C=[2*zeta1,0;0,2*zeta2*mu]+epsilon*[0 c12;c21,0];
    % eig(K,M)
    %
    Preload=[0;0;omega;0];
    %
    A_matrix=[zeros(2,2),eye(2);-MHH\K,-MHH\C];
    % eig(A_matrix)
    % sweep the list
    
    % time integration
    [flag,tout,yout,yeout0,teout1,yeout1,ieout]=DirectInt(m_c,res,A_matrix,...
        InitCond,tspan,MHH,Preload);
%     figure(1)
%     plot(tout,yout(:,1))
%     figure(2)
%     plot(tout,yout(:,2))
%     figure(3)
%     plot(tout,yout(:,3))
%     figure(4)
%     plot(tout,yout(:,4))
    %
    t_kept=T-100;
    abstract_ind=(teout1>t_kept);
    if any(abstract_ind)
        try
%         yeout0(end,:)=[];
        p1_section_in=yeout0(abstract_ind,3);
        p1_section_out=yeout1(abstract_ind,3);
        p2_section_in=yeout0(abstract_ind,4);
        p2_section_out=yeout1(abstract_ind,4);
        
        m1_v_sec=[unique(p1_section_in);unique(p1_section_out)];
        m2_v_sec=[unique(p2_section_in);unique(p2_section_out)];
        catch ME
            disp('Wrong message')
            keyboard
        end
    else
        abstract_ind=(tout>0.9*(tout(end)));
        m1_v_sec=[max(unique(yout(abstract_ind,3)));min(unique(yout(abstract_ind,3)))];
        m2_v_sec=[max(unique(yout(abstract_ind,4)));min(unique(yout(abstract_ind,4)))];
    end
    % get rid ofduplicate
    m1_temp=sort(m1_v_sec,'descend')';
    front_serial=m1_temp(1:end-1);
    latter_serial=m1_temp(2:end);
    adjacent_dis=abs(front_serial-latter_serial);
    m1_temp=[m1_temp(1),latter_serial(adjacent_dis>1e-6)];
    m2_temp=sort(m2_v_sec,'descend')';
    front_serial=m2_temp(1:end-1);
    latter_serial=m2_temp(2:end);
    adjacent_dis=abs(front_serial-latter_serial);
    m2_temp=[m2_temp(1),latter_serial(adjacent_dis>1e-6)];
    %
    m1_v_temp=[bif_para,m1_temp];
    m2_v_temp=[bif_para,m2_temp];
    %
    length_m1_v(j)=length(m1_v_temp);
    length_m2_v(j)=length(m2_v_temp);
    %
    m1_v=[m1_v,m1_v_temp];
    m2_v=[m2_v,m2_v_temp];
%     catch ME 
%     end
end
toc
% post processing
n_col=max(length_m1_v);
m1_diagram=zeros(length(length_m1_v),n_col);
temp_indx=1;
for j=1:length(length_m1_v)
    m1_diagram(j,1:length_m1_v(j))=m1_v(temp_indx:temp_indx+length_m1_v(j)-1);
    m1_diagram(j,length_m1_v(j)+1:end)=nan;
    temp_indx=temp_indx+length_m1_v(j);
end
%
n_col=max(length_m2_v);
m2_diagram=zeros(length(length_m2_v),n_col);
temp_indx=1;
for j=1:length(length_m2_v)
    m2_diagram(j,1:length_m2_v(j))=m2_v(temp_indx:temp_indx+length_m2_v(j)-1);
    m2_diagram(j,length_m2_v(j)+1:end)=nan;
    temp_indx=temp_indx+length_m2_v(j);
end
filename=sprintf('%dbifurcation_S.mat',yyyymmdd(datetime));
save(filename,'m1_diagram','m2_diagram')



%

figure
plot(m1_diagram(:,1),m1_diagram(:,2:end),'b.','markersize',8)
xlabel(x_b,'interpreter','latex')
ylabel('$\dot{x_{1}}$','interpreter','latex')
set(gca,'fontname','times new roman')
xlim([0 0.03])
figure
plot(m2_diagram(:,1),m2_diagram(:,2:end),'b.','markersize',8)
xlabel(x_b,'interpreter','latex')
ylabel('$\dot{x_{2}}$','interpreter','latex')
set(gca,'fontname','times new roman')
xlim([0 0.03])





















% sub_functions
function [flag,tout,yout,yeout0,teout1,yeout1,ieout]=DirectInt(m_c,res,A_matrix,InitCond,tspan,MHH,Preload)
flag=0;
fs=100;
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
f_sk= @(t,y,A_matrix,Preload) f_ls(t,y,A_matrix,Preload)+ (hx_(t,y)*f_ls(t,y,A_matrix,Preload))*Wx_(t,y)/(1+res);

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
%         f_ls(t,y0,A_matrix,Preload)
%         f_sk(t,y0,A_matrix,Preload)
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


