function [x] = dichotomy(F, a, b ,Tol, debug)
%> get the funciton's root
N = ceil( log2( (b-a)/Tol ) );
fa = F(a);
fb = F(b);
x_list = [a, b];
y_list = [fa, fb];
if (fa*fb>0)
    disp('The root is not in this range!')
    return
else
    if debug    
        FIG = figure;
%         haxes = axes;
    end
    
    n   = 1;
    while n < N +2  && abs(b-a)>Tol
        x = (a + b)/2;
%         x = a - (b-a)/(fb-fa)*fa;
        
        fa = F(a);
       
        fx = F(x);
         if debug
             close all
%              FIG = figure;
             figure
             plot( [a b], [fa fb], 'ro')
             hold on
             plot(x_list(1:2), [0 0], 'b-')
             plot(x,fx, 'k*')
             x_list = [x_list, x];
             y_list = [y_list, fx];
         end
        %
        if fa*fx>0
            a = x;
        else
            b = x;
        end
        n = n + 1;
    end
end

end