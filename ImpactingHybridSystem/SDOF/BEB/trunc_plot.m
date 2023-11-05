% 
function [ytop,ybot,xright,CASE,BEB2D]=trunc_plot(tout,yout,teout,trunc,CASE,SIGNCASE,fig_num,BEB2D)

public = BEB2D.public;
fs =SIGNCASE.fs;
% get the event mark location
index_1 = abs(tout'-teout)< 0.01/fs;
% the result may be a binary vector or a matrix
if size(index_1,1) ==1 %case of a row vector
    index1 = index_1;
else %case of a binary matrix
   index1 =sum(index_1)';
end
% then finnaly convert it to a collumn vector
if size(index1,1)==1
    index1 =index1';
end

% according to given time T to trunct the last time history
trun_sequence = zeros(1,length(tout))';
trun_sequence(floor(end-fs*trunc):end)=1;
trunc_loc = trun_sequence>0;
% 
yout = yout(trunc_loc,:);
trun_x = yout(:,1);
trun_y = yout(:,2);
%
xright = max(trun_x);
ytop =max(trun_y);
ybot = min(trun_y);




% update the event location in new series
index1 = index1 & trun_sequence;
% count the discarded points
n_q = length(tout)- fs*trunc-1;
loc= find(index1>0);
% update the new location 
loc = loc - n_q;

% count number of remaining hits
n_h = length(loc);
%
% plot in the appointed figure number
SIGNCASE.segnum = n_h-2;
State =[trun_x(1:loc(1)),trun_y(1:loc(1))]';
SIGNCASE.seg00.X = trun_x(1:loc(1));
SIGNCASE.seg00.Y = trun_y(1:loc(1));
SIGNCASE.seg00.S = State;
for  i=2:n_h-1
    segname = strcat('seg',num2str(i-1));
    
    State =[trun_x(loc(i-1)+1:loc(i)),trun_y(loc(i-1)+1:loc(i))]';
    
    eval(['SIGNCASE.', segname, '.X = trun_x(loc(i-1)+1:loc(i));'])
    eval(['SIGNCASE.', segname, '.Y = trun_y(loc(i-1)+1:loc(i));'])
    
    eval(['SIGNCASE.', segname, '.S = State;'])
    
    %
end
% the last outcoming seg
State =[trun_x(loc(n_h)+1:end),trun_y(loc(n_h)+1:end)]';
SIGNCASE.seg01.X = trun_x(loc(n_h)+1:end);
SIGNCASE.seg01.Y = trun_y(loc(n_h)+1:end);
SIGNCASE.seg01.S = State;


% virtual trajectory implement by interpolation 
% unfinished 
% virtual segment num is n_h -2
% for  j =1 : n_h-2
    % first 
% end

%----------------------------------------------------------------------%
% inherit from the case
if sign(SIGNCASE.mu)>0
    CASE.POS = SIGNCASE;
else
    CASE.MINUS = SIGNCASE;
end

% PLOT : change the time history via the observing function
figure(fig_num)
S_temp = SIGNCASE.seg00.S;
mu = SIGNCASE.mu;
[X,Y] = Observing_filter(public,CASE,S_temp,mu);
plot(X,Y,'-','color',public.plot.addmissible_ob.color,'linewidth',public.plot.addmissible_ob.linewidth)
hold on
for  i=2:n_h-1
    segname = strcat('seg',num2str(i-1));
    eval(['S_temp =','SIGNCASE.', segname, '.S;'])
    [X,Y] = Observing_filter(public,CASE,S_temp,mu);
    plot(X,Y,'-','color',public.plot.addmissible_ob.color,'linewidth',public.plot.addmissible_ob.linewidth)
end
S_temp = SIGNCASE.seg01.S;
[X,Y] = Observing_filter(public,CASE,S_temp,mu);
plot(X,Y,'-','color',public.plot.addmissible_ob.color,'linewidth',public.plot.addmissible_ob.linewidth)

%

eval(['BEB2D.',strcat('CASE',num2str(CASE.num)),'=CASE';]);
end