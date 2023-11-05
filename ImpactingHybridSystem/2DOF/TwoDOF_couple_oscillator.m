% 2DOF Impact Oscillators
clc
clear
close all
% Parameter Definition

%
tspan=0:0.02:4*pi;
zeta1=-0.1; zeta2=0.1;
mu=2;   % w2/w1
eta=5;   % m2/m1
epsilon=0.1; % the perturbation

m21=1;m12=m21*eta;
c21=-1;c12=c21*eta;
k21=-1;k12=k21*eta;
X1_0=0;dX1_0=1;X2_0=0;dX2_0=1;

% derived parameters
a1=-zeta1;a2=-zeta2*mu;b1=sqrt(1-zeta1^2);b2=mu*sqrt(1-zeta2^2);
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
Am11=Amplitude(a2-a1,b2+b1,-A1,A2);
Am12=Amplitude(a2-a1,b2+b1,A2,A1);
Am13=Amplitude(a2-a1,b2-b1,-A1,A2);
Am14=Amplitude(a2-a1,b2-b1,-A2,-A1);
Am15=Amplitude(a2-a1,b2+b1,A1,-A2);
Am16=Amplitude(a2-a1,b2+b1,A2,A1);
Am17=Amplitude(a2-a1,b2-b1,A1,-A2);
Am18=Amplitude(a2-a1,b2-b1,A2,A1);

if X2_0~=0
    Am20=X2_0/cos(phi_2);
else
    Am20=dX2_0/(a2*cos(phi_2)-b2*sin(phi_2));
end
Am21=Amplitude(a2-a1,b2+b1,-B1,-B2);
Am22=Amplitude(a2-a1,b2+b1,-B2,B1);
Am23=Amplitude(a2-a1,b2-b1,-B1,B2);
Am24=Amplitude(a2-a1,b2-b1,-B2,-B1);
Am25=Amplitude(a2-a1,b2+b1,B1,B2);
Am26=Amplitude(a2-a1,b2+b1,-B2,B1);
Am27=Amplitude(a2-a1,b2-b1,B1,-B2);
Am28=Amplitude(a2-a1,b2-b1,B2,B1);

% time history function
X10=@(t) Am10*exp(a1.*t).*cos(b1*t+phi_1);
X20=@(t) Am20*exp(a2.*t).*cos(b2*t+phi_2);

% X11=@(t) exp(a2*t).*(Am11.*(cos(b2*t+phi_1+phi_2)-cos(phi_1+phi_2))+...
%     Am12.*(sin(b2*t+phi_1+phi_2)-sin(phi_1+phi_2))+...
%     Am13.*(cos(b2*t+phi_2-phi_1)-cos(phi_2-phi_1))+...
%     Am14.*(cos(b2*t+phi_2-phi_1)-cos(phi_2-phi_1)));
% 
% X21=@(t) exp(a1*t).*(Am21.*(cos(b1*t+phi_1+phi_2)-cos(phi_1+phi_2))+...
%     Am22.*(sin(b1*t+phi_1+phi_2)-sin(phi_1+phi_2))+...
%     Am23.*(cos(b1*t+phi_1-phi_2)-cos(phi_1-phi_2))+...
%     Am24.*(cos(b1*t+phi_1-phi_2)-cos(phi_1-phi_2)));

X11=@(t) 0.5/b1 * (...
    exp(a1*t).*(Am11*cos(b1*t+phi_1-phi_2)+Am12*sin(b1*t+phi_1-phi_2)+...
    Am13*cos(b1*t+phi_1+phi_2)+Am14*sin(b1*t+phi_1+phi_2))+...
    exp(a2*t).*(Am15*cos(b2*t-phi_1+phi_2)+Am16*sin(b2*t-phi_1+phi_2)+...
    Am17*cos(b2*t+phi_1+phi_2)+Am18*sin(b2*t+phi_1+phi_2))...
                   );

X21=@(t) 0.5/b2 * (...
    exp(a1*t).*(Am21*cos(b1*t+phi_1-phi_2)+Am22*sin(b1*t+phi_1-phi_2)+...
    Am23*cos(b1*t+phi_1+phi_2)+Am24*sin(b1*t+phi_1+phi_2))+...
    exp(a2*t).*(Am25*cos(b2*t-phi_1+phi_2)+Am26*sin(b2*t-phi_1+phi_2)+...
    Am27*cos(b2*t+phi_1+phi_2)+Am28*sin(b2*t+phi_1+phi_2))...
                   );
X11(0)
X21(0)
% Original System


M=[1,0;0,1]+epsilon*[0 m12;m21,0];
K=[1,0;0,mu^2]+epsilon*[0 k12;k21,0];
C=[2*zeta1,0;0,2*zeta2*mu]+epsilon*[0 c12;c21,0];
% eig(K,M)
%
omega=0;
%
A_matrix=[zeros(2,2),eye(2);-M\K,-M\C];
% eig(A_matrix)

F=@(t,y,A_matrix,omega) A_matrix*y+[0;0;omega;0];

y0=[X1_0;X2_0;dX1_0;dX2_0];
[t,y]=ode45(@(t,y)F(t,y,A_matrix,omega),tspan,y0);



%
figure
plot(t,y(:,1),'k--','linewidth',2)
hold on
plot(t,X10(t),'r','linewidth',1.2)
legend('Mass1-Num','Mass2-Ana')


figure
plot(t,y(:,2),'k--','linewidth',2)
hold on
plot(t,X20(t),'r','linewidth',1.2)
legend('Mass2-num','Mass2-Ana')

% cor1=X11(t);
% cor2=X21(t);

figure
plot(t,X11(t),'k','linewidth',1.2)
hold on
plot(t,X21(t),'r','linewidth',1.2)
legend('Per_Mass1','Per_Mass2')


% compare with the perturbation method and the direct num method

figure
subplot(121)
plot(t,y(:,1),'k--','linewidth',2)
hold on
plot(t,X10(t)+epsilon*X11(t),'r-','linewidth',1.2)
legend('Mass1-DN','Mass1-Pert')
%
subplot(122)
plot(t,y(:,2),'k--','linewidth',2)
hold on
plot(t,X20(t)+epsilon*X21(t),'r-','linewidth',1.2)
legend('Mass2-DN','Mass2-Pert')

function Am=Amplitude(A,B,C,D)

%
Am=(C*A+D*B)/(A^2+B^2)/2;

end

