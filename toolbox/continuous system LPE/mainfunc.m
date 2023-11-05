clear;
clc;
% 连续系统的LPE计算
SIGMA_list =  0.5:0.001:0.9;
R = 1; 
BETA = 1;

 %       
dt=0.01;
time=500;
y0=[0.01;0.01;0.01];
tspan=0:dt:time;
cut=20000;
k=1;
err=1e-3;
figure(1)
LP=zeros(numel(SIGMA_list),1);
for i=1:numel(SIGMA_list)
    SIGMA=SIGMA_list(i); 
    params=[SIGMA,R,BETA];
    [t,y]=ode45(@(t,y)dydt(t,y,params),tspan,y0);
    Y=y;
    clf(1)
    plot(Y(cut:end,1),Y(cut:end,2))
    Y1=Y(cut:end,1);
    Y2=Y(cut:end,2);
    Y3=Y(cut:end,3);
    loc=find(abs(diff(Y1)./diff(Y2)-k)<=err&diff(Y2)>=0);
    hold on
    if isempty(loc)~=1
        scatter(Y1(loc+1),Y2(loc+1));hold off
        store.Y1{i}=Y1(loc+1);
        store.Y2{i}=Y2(loc+1);
        store.Y3{i}=Y3(loc+1);
    else
        store.Y1{i}=nan;
        store.Y2{i}=nan;
        store.Y3{i}=nan;
    end
    pause(0.1)
    % get LPE
    [T,Res]=lyapunov(3,@Dydt,@ode45,0,0.5,500,y0,0,params);
    for k=1:3
    LP(i,k)=sum(Res(end-200:end,k))/201;
    end
end
figure(2)
hold on
for i=1:numel(SIGMA_list)
    plot(SIGMA_list(i),store.Y2{i},'Color','b','Marker','.')
end

figure(3)
hold on
for i=1:numel(SIGMA_list)
    hold on
    plot(SIGMA_list,LP(:,1),'Color','r','Marker','.')
    plot(SIGMA_list,LP(:,2),'Color','y','Marker','.')
    plot(SIGMA_list,LP(:,3),'Color','g','Marker','.')
end

%_______________________________________________________%

function dydt=dydt(t,X,params)
SIGMA = params(1); R = params(2); BETA = params(3);
x=X(1); y=X(2); z=X(3);
dydt=    [y;...
            z;...
            -SIGMA*z-y-R*x^2-BETA*x];
end
function dydt=Dydt(t,X,params)
SIGMA = params(1); R = params(2); BETA = params(3);
x=X(1); y=X(2); z=X(3);

Q= [X(4), X(7), X(10);
    X(5), X(8), X(11);
    X(6), X(9), X(12)];
dydt=zeros(9,1);

dydt(1:3)=    [y;...
            z;...
            -SIGMA*z-y-R*x^2-BETA*x];
Jac=[0,1,0; 0,0,1; -2*R*x-BETA, -1,-SIGMA];

dydt(4:12)=Jac*Q;
        
end
