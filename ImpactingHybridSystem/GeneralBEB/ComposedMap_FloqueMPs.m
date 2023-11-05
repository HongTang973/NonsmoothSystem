function [MainV,LPE]=ComposedMap_FloqueMPs(x_p,T,fs,F_ls,RMap,EventFun,PlotInd)
% For a general Hybrid system with a flow region and the reset map on the
% discontinuity surface, if there is a LCO with
% x_0: incoming point
% x_00: outcoming point
% x_p : point on the poincare section
% T: the period of this orbit
% F_ls: function handle of the flow
% RMap: the reset map when the flow hit the boundary
% EventFun: the function to detect the event of hitting the boundary
% PlotInd: 1x5 vector, with the first element  specifying the obsevor 
% next two variables to plot out phase, last integers to specify the
% figure number to open
% *********************************************************%
%---------------------------AIM----------------------------%
% Calculate the Floque Multipliers around the LCO to check %
% the stability of the LCO                                 %
% *********************************************************%

% n_s : the dimension of the system
n_s = length(x_p);
% scaler to make the change very tiny, close to limit of differentiation
Scale = 1000000;
delta=rand(n_s,n_s);
delta_xp = delta/Scale;
% define the num of loop and store the LPE
N_loop=100;
LPE=[];
for  k = 1:N_loop
    Deviation0=[];
    % define the different initial condition with 
    % different small pertuebations 
    for  i=1:n_s; eval(['x',num2str(i), '= x_p+delta_xp(:,',num2str(i),');']);  end
    % start integrating with different initial conditions
    for  j=1:n_s 
        eval( ['X',num2str(j),'=composedMAP(x',num2str(j),',F_ls,RMap,EventFun,T,fs,PlotInd);']);
        %  get the deviation vectors regarding to the initial point
        Deviation0=[Deviation0,eval(['X',num2str(j),'-x_p'])];
    end
   % get the LPE and Construct new orthonormal basis by gram-schmidt
   % orthogonalization 
    [lpe,delta_xp]= LPE_GSC(Deviation0,delta_xp);
    LPE=[LPE,lpe];
    % Get the new run of flow  with new initial conditions
end
% 
MainV = delta_xp;
for k=1:n_s
        % to make the initial purtbation vectors very tiny
         MainV(:,k)=MainV(:,k)/norm(MainV(:,k));
end

% -------- BUILT IN FUNCTIONS AFTER THE MAIN FUNCTION-----------%
%  integration function
function  [X]=composedMAP(y0,F_ls,RMap,EventFun,T,fs,PlotInd)
% close all 
refine=1;
EventFun(0,y0)
options = odeset('Events',@(t,y) EventFun(t,y),...
    'RelTol',1e-12,'AbsTol',1e-12,'Refine',refine);
%
timespan=[0:1/fs:T];
% Solve until the first terminal event: hit the boundary.
[t1,y1,te,ye,~] = ode45(@(t,y) F_ls(t,y),timespan,y0,options);

% reset map
y10=RMap(ye');
% restart integrating until the T.
timespan=[te:1/fs:T];
[t2,y2,~,~,~] = ode45(@(t,y) F_ls(t,y),timespan,y10,options);
% Get the ending point after the evolution time T
X=y2(end,:)';

% plot the process about the time history and the phase
figure(PlotInd(4))
plot(t1,y1(:,PlotInd(1)))
hold on
plot(t2,y2(:,PlotInd(1)))
D=[y1(:,PlotInd(2));y2(:,PlotInd(2))];
V=[y1(:,PlotInd(3));y2(:,PlotInd(3))];
%
figure(PlotInd(5))
plot(D,V)
hold on
plot(y0(PlotInd(2)),y0(PlotInd(3)),'r*')
plot(X(PlotInd(2)),X(PlotInd(3)),'bo')


%
function [lpe,New_delta_xp]= LPE_GSC(Deviation0,delta_xp)
% Construct new orthonormal basis by gram-schmidt orthogonalization
N= size(Deviation0,1);
lpe=zeros(N,1);
%
scale=1000000;
% project the Derivation vector on the initial purtbation vectors
for i = 1:N; lpe(i)=(Deviation0(:,i)'*delta_xp(:,i))/norm(delta_xp(:,i))^2; end
%
New_delta_xp = Deviation0;
New_delta_xp(:,1)= New_delta_xp(:,1)/norm(New_delta_xp(:,1));
    for j=2:N
        for k=1:j-1 
            %get the extent of containing the previous vectors
            gsc= Deviation0(:,j)'*New_delta_xp(:,k)/norm(New_delta_xp(:,k))^2;  
            % subtruct the components
            New_delta_xp(:,j)=New_delta_xp(:,j)-gsc*New_delta_xp(:,k);
        end        
%         Basis(:,j)'*Basis(:,1)
    end
    for k=1:N
        % to make the initial purtbation vectors very tiny
         New_delta_xp(:,k)=New_delta_xp(:,k)/norm(New_delta_xp(:,k))/scale;
    end

%------------END OF THE MAIN FUNCTION------%
