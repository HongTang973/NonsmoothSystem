% 连续系统的LPE计算
SIGMA = 10; R = 28; BETA = 8/3;
params=[SIGMA,R,BETA];
ystart=[0 1 0]; % 不同的初始条件可能有不同的吸引子 不同的吸引子的LPE也不一样
[T,Res]=lyapunov(3,@lorenz_ext,@ode45,0,0.5,200,ystart,10,params); 
figure
plot(T,Res)