function [roots, fvals] =get_roots_Det(A,R,C, tspan,debug)
D1          =eig(A);
Omega       = max(abs(D1));

fs          = 1e4;
if log10(tspan(2) - tspan(1))<0
    fs = fs /10^floor(log10(tspan(2) - tspan(1)));
end
delta       = 1/fs;
T           = tspan(1):delta:tspan(2);
%
Det         = @(t) det(R*expm(A*t) - eye(length(C)));

%> initialize the container of roots
roots_ind   = [];
%> set the default sign as positive
temp_sign   = 1;
fval_list   = zeros(1, length(T));
%> calculate the value of the det function and detect the sign change
for i=1:length(T)
    F_1= Det(T(i));
    if temp_sign*sign(F_1) <=0  %> detect the sign change
        roots_ind = [roots_ind, i];
    end
    %> store the value of current time step and its sign
    fval_list(i) = F_1;
    temp_sign = sign(F_1);
end

%
if debug
    FIG =figure;
    plot(T,fval_list)
    hold on
    plot(tspan,[0 0], 'r--')
    xlabel('T')
    ylabel('val')
    adj_plot_theme_I(FIG)
end
%
roots = [];
fvals = [];

if isempty(roots_ind)
    disp('No roots in the specified range!')
else
    for j = 1: length(roots_ind)
        %> T(roots_ind)
%         
        try 
            temp_root = dichotomy( Det, T( roots_ind(j) -1 ),  T( roots_ind(j) ), 1e-9, 0);
        catch ME
            temp_root = fzero( Det, T( roots_ind(j) ) );
        end
        %> use the dichotomy method to get the roots
        
        roots = [roots, temp_root];
        fvals = [fvals,Det(temp_root)];
     
    end
end
%> get rid of the fake roots
        fvals(roots<1e-4) =[];
        roots(roots<1e-4) =[];
%> get rid of the duplicate roots
if ~isempty(roots)
    [roots, index] =sort(roots,'ascend');
     fvals = fvals(index);
    front_serial =roots(1:end-1);
    latter_serial=roots(2:end);
    latter_serial_f = fvals(2:end);
    adjacent_dis=abs(front_serial-latter_serial);
    roots =[roots(1),latter_serial(adjacent_dis>1e-6)];
    fvals = [fvals(1), latter_serial_f(adjacent_dis>1e-6)];
end

end