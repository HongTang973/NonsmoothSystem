% Detecting the existence of LCO when boundary equilibrium transferring to
% virtual Equilibrium for general n dimension impacting hybrid system 
% limitations: single impact; single discontinuity surface

% KEY INPUTs:
% A : the linear region system's matrix 
% R : the reset map matrix, R(x) = R*x;
% C : the selecting vector to define the discontinuity surface C*x=0
% 

% TRY with a Valve_QWM  problem: model dimension - 5
clc
clear
close all
% if equi_type == 'admissible'
%         pre = -1;
%  elseif equi_type == 'pseudo'
%         pre = 1;
%  end
equi_type = 1;

% define the linear part
A =[ 0    1.0000         0         0         0;
   37.7917   -0.0500    2.7298    2.7298         0;
   -0.0249         0   -0.0011   -0.0011   -0.0001;
    0.0352         0    0.0016    0.0016    0.4740;
   -0.1290 -280.8643   -0.0058 -221.4821   -5.9348];


% define the reset map
pr=0.78; %coef of restitution
C  = [1 0 0 0 0];
% vx_= @(t,y)  C*A*(y-state0); 
B= -(1+pr)*[0;  1; 0 ; 0; 0];
% define initial condition: grid on the positive set
MW = B*C*A;
R = eye(length(C))+MW;
% [VV,DD] =eig(R)
% C*A*VV
%  OP=P*R*inv(P)
% GG=R'*R -eye(3)
% [v,d]=eig(GG)
% eig(GG(2:3,2:3))
% %
A_s = (eye(5)-(B*C*A)/(C*A*B))*A;
eig(A_s)
[V1,D1]=eig(A);
EA = @(T) real(V1*diag(exp(diag(D1)*T))*inv(V1));
sign_V = @(Y) sign(C*A*Y);

% determine the sampling frequency 
Omega = max(abs(diag(D1)));
fs = 10*ceil(2*(Omega/2/pi));



% for different Evolution time
a =0;
b = 2;
delta = 0.001;
T = a*pi:delta:b*pi;
%
MAX = zeros(1,length(T));
F_1 = zeros(1,length(T));
V_sign = zeros(1,length(T));
LOCI= zeros(size(A,1),length(T));
vector= zeros(size(A,1),length(T));
for i=1:length(T)
%     [V_sign(i),LOCI(:,i),MAX(i),F_1(i)] = solve_classification(T(i),EA,sign_V);

    [V_sign(i),LOCI(:,i),MAX(i),vector(:,i),F_1(i),COND] = LCO_Det_search(T(i),R,A,C,equi_type);

%     [V_sign(i),LOCI(:,i),MAX(i),vector(:,i),F_1(i)] = LCO_detecting_line_search(T(i),R,EA,sign_V,C,equi_type);
end

%% make the plot
velocity = C*A*vector;
figure(1)
plot( T/pi,F_1,'r.','linewidth',1.2,'displayname','F1')
hold on
% plot( T/pi,V_sign,'-','linewidth',1.2,'displayname','Sign of Velocity')
plot([a b],[0 0],'b-','displayname','zero line')
% plot( T/pi,velocity,'-','linewidth',1.2,'displayname','Sign of Velocity')
% plot([a b],[0 0],'b-','displayname','zero line')
% plot( T/pi,MAX,'k-','linewidth',1.0,'displayname','MAX-1')
plot( T/pi,LOCI(1,:),'-','linewidth',1.0,'displayname','E-1')
plot( T/pi,LOCI(2,:),'-','linewidth',1.0,'displayname','E-2')
plot( T/pi,LOCI(3,:),'-','linewidth',1.0,'displayname','E-3')
legend
ylim([-40 40])
xlabel('\tau /\pi')
ylabel('value')
set(gca,'fontname','Times New Roman')

%%
B0 =B;
C*A*B0
[A_c,B_c,PZ] = companion_form(A,-B0,C')
% eig(A_c)
% B_c = B((PZ\P_l3)/norm(PZ\P_l3),(PZ\V_l3)/norm(PZ\V_l3))
% stability of PE
A_s_c = (eye(5)-(B_c*C*A_c)/(C*A_c*B_c))*A_c;
eig(A_s_c)
A_s = (eye(5)-(B0*C*A)/(C*A*B0))*A;
eig(A_s)
pec = pinv([C;C*A_c;A_s_c])*[C';0;0]
pe = pinv([C;C*A;A_s])*[C';0;0]

PZ*pec
C*A_c*A_c*B_c
%% post processing

% there is a/some cross points
index0 =sign(F_1);
index1 = abs(diff(index0))>0;
% filter the singularity case
index_s = abs(diff(F_1))/delta < (1/delta);
index1 =index1 & index_s;
index2 = [index1,0];
index3 = [0,index1];
index2=find(index2==1);
index3=find(index3==1);
% section partia
ratio = abs(F_1(index2))./(abs(F_1(index2))+abs(F_1(index3)));
T_chosen =(1-ratio).*T(index2)+ratio.*T(index3);
MAX_chosen = (1-ratio).*MAX(index2)+ratio.*MAX(index3);
F1_chosen = (1-ratio).*F_1(index2)+ratio.*F_1(index3);
Sign_chosen =(1-ratio).*V_sign(index2)+ratio.*V_sign(index3);
%          figure(1)
% plot(T_chosen/pi,Sign_chosen,'bo')
%     %
% index4 = (Sign_chosen == 1) & (abs(F1_chosen)<1e-3) & (abs(MAX_chosen)<1e-3);
% T_chosen  =T_chosen(index4);
% F1_chosen =F1_chosen(index4);


%

if ~isempty(T_chosen)
    LCO = zeros(length(C),length(T_chosen));
    if length(T_chosen)>2
    left = T_chosen(1)-0.2*sum(abs(T_chosen(end)-T_chosen(1)));
    right =T_chosen(end)+0.2*sum(abs(T_chosen(end)-T_chosen(1)));
    else
    left = T_chosen(1)-0.1*sum(abs(T_chosen));
    right =T_chosen(end)+0.1*sum(abs(T_chosen));
    end
    for i=1:length(T_chosen)
        % get the LCO with statevariables
        [~,~,~,LCO(:,i),~] = LCO_detecting_line_search(T_chosen(i),R,EA,sign_V,C,equi_type);
        figure(1)
        plot(T_chosen(i)/pi,F1_chosen(i),'bo','displayname',['the ',num2str(i),'th candidate'])
    end
    figure(1)
    xlim([left right]/pi)
    legend('location','best')
    disp([num2str(length(T_chosen)),' LCO(s) found!'])
    LCO
    disp('with period:')
    T_chosen/pi
    for i=1:length(T_chosen)
        %check the floque multipliers of the foud LCO
        [Mono_p,Salt_p]=Floque_Multipliers(T_chosen(i),EA,LCO(:,i),R,A,C)         
        if max(abs(Salt_p))>1+1e-6
            disp([' LCO',num2str(i),' is unstable!'])
        else
            disp([' LCO',num2str(i),' is stable!'])
        end
    end
else
    disp('No LCO found!')
end
keyboard
%
if ~isempty(T_chosen)
    tspan = [0 4];
    for i =1:length(T_chosen)
    InitCond = LCO(:,i);
    [tout,yout,yeout0,teout,yeout,ieout]=...
        Single_DS_Impacting_Hybrid_system_integration(A,R,C,InitCond,tspan,fs,equi_type);
    eval(['tout',num2str(i),'=tout;']);
    eval(['yout',num2str(i),'=yout;']);
    eval(['yeout0',num2str(i),'=yeout0;']);
    eval(['teout',num2str(i),'=teout;']);
    eval(['yeout',num2str(i),'=yeout;']);
    eval(['ieout',num2str(i),'=ieout;']);
    eval(['delta_t',num2str(i),'=teout',num2str(i),'(2:end)-teout',num2str(i),'(1:end-1);']);
    %Period
    eval(['Period',num2str(i),'=delta_t',num2str(i),'/pi;']);
    figure
    plot(tout,yout(:,C>0),'r-','linewidth',1.4)
    title(['LCO ',num2str(i),'''s stability'])
    set(gca,'fontname','times new roman','fontsize',12)
    xlabel('t/s')
    figure
    plot3(yout(:,1),yout(:,2),yout(:,3),'k-','linewidth',1.2)
    title(['LCO ',num2str(i),'''s phase portrait'])
    set(gca,'fontname','times new roman','fontsize',12)
    grid on
    end
end
% tspan = [0 2000];
% pert = [-0.121018860540516;0.642369009705563;0.756780345122262];
% pert(1)=0;
% epsilon=0.00;
% [tout,yout,yeout0,teout,yeout,ieout]=...
%         Single_DS_Impacting_Hybrid_system_integration(A,R,C,LCO(:,1)+epsilon*pert,tspan,fs);
%     figure
% plot(tout,yout(:,1))

[tout,yout,yeout0,teout,yeout,ieout]=...
        Single_DS_Impacting_Hybrid_system_integration(A,R,C,-LCO(:,1),[0 100],fs,-equi_type);
    figure
    plot(tout,yout(:,C>0),'r-','linewidth',1.4)