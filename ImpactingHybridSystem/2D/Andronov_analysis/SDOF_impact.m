% 1DOF function
% nondimensionalized
clc
clear
close all
% mu=0.1;
% r=1/0.7;
mu=-0.1;
r=0.7;

% Innitial condition
x_0=-1;
y_0=0.1;
%
x=@(x_0,y_0,tau)  exp(-mu.*tau).*(x_0.*cos(tau)+(y_0.*sqrt(1+mu^2)+mu.*x_0).*sin(tau));
y=@(x_0,y_0,tau)  exp(-mu.*tau).*(y_0.*cos(tau)-(x_0.*sqrt(1+mu^2)+mu.*y_0).*sin(tau));

F=@(tau) 1-exp(mu*tau)*(cos(tau)-mu*sin(tau));
phi=@(tau,mu) 1-exp(mu*tau)*(cos(tau)-mu*sin(tau));

% set preparation for root finding
    options = optimoptions('fsolve','Display','off','FunctionTolerance',1e-9);
% get the critical V
t0=1.5*pi;
[dt1,fval] = fsolve(F,t0,options);
dt1/pi;
V=-exp(mu*dt1)*phi(dt1,-mu)/sqrt(1+mu^2)/sin(dt1)
% y_0=8.35775*V

% y_0=1
V_=exp(-mu*dt1)*phi(dt1,mu)/sqrt(1+mu^2)/sin(dt1)
% time history
 T=1000;
dt=0.01;
yout=[];
xout=[];
tout=[];
num=10;
t_final=0;
while t_final<T
        dt_temp=[];
        X=@(tau)  exp(-mu.*tau).*(x_0.*cos(tau)+(y_0.*sqrt(1+mu^2)+mu.*x_0).*sin(tau))+1;
%         t_net=0:0.01:2*pi;
%         x_net=X(t_net);
%         figure(1)
%         plot(t_net,x_net,'k-','linewidth',1.2)
        for i=1:num % search roots among grid seeds 
            t0=i/num*2*pi;
            [dt1,fval] = fsolve(X,t0,options);
%             [dt2,~]=instant_solve(beta,dot_beta,mcoef,eigen_value,t0);
            dt_temp=[dt_temp;real(dt1),fval];
        end
        
        % erase the trival roots
        dt_temp(dt_temp<1e-3)=[];
        dt_temp(dt_temp>6.18)=[];
        dt_temp=(sort(real(dt_temp)));
        if isempty(dt_temp)||abs(X(dt_temp(1)))>1e-3
            t_list=0:0.1:T;
            xout1=x(x_0,y_0,t_list);
            yout1=y(x_0,y_0,t_list);
            tout1=t_list+t_final;
            xout=[xout,xout1];
            yout=[yout,yout1];
            tout=[tout,tout1];
            break
        else
            % hit the boundary
            t1=dt_temp(1);
            t_list=0:0.1:t1;
            t_list=unique([t_list,t1]);
            xout1=x(x_0,y_0,t_list);
            yout1=y(x_0,y_0,t_list);
            tout1=t_list+t_final;
            xout=[xout,xout1];
            yout=[yout,yout1];
            tout=[tout,tout1];
%             xout(end)
            %
            
            y_0=-r*yout(end);
            x_0=-1;
            t_final=t_final+t1;
        end
        
end


figure
plot(tout,xout,'k-','linewidth',1.2)
set(gca,'fontname','times new roman','fontsize',12)
figure
plot(tout,yout,'k-','linewidth',1.2)
set(gca,'fontname','times new roman','fontsize',12)
figure
plot(xout,yout,'k-','linewidth',1.2)
set(gca,'fontname','times new roman','fontsize',12)


