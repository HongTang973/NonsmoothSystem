% define the input
% t0
% td
% c
% rhs = 2*A_t/A_s
t0 =100;
td =5;
c = 1;
rhs =0.1;

% define a function 
F = @(T) exp(-(T-td).^2/2/t0)- c*exp(-T.^2/2/t0) -rhs;
%
delta =0.1;
sample_T = 0:delta:10;
v_F  =  F(sample_T);

% there is a/some cross points
index0 =sign(v_F);
index1 = abs(diff(index0))>0;
% filter the singularity case
index_s = abs(diff(v_F))/delta < (1/delta);
index1 =index1 & index_s;
index2 = [index1,0];
index3 = [0,index1];
index2=find(index2==1);
index3=find(index3==1);
% section partia
ratio = abs(v_F(index2))./(abs(v_F(index2))+abs(v_F(index3)));
T_chosen =(1-ratio).*sample_T(index2)+ratio.*sample_T(index3);
F1_chosen = (1-ratio).*v_F(index2)+ratio.*v_F(index3);
if ~isempty(T_chosen)

figure
plot(sample_T,v_F,'displayname',['the function line'])
hold on
plot(T_chosen,F1_chosen,'bo','displayname',['the solution point'])
legend('location','best')
disp(['In this case, your solution is :'])
T_chosen
end