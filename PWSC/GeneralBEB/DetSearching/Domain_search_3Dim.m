% using 2D searching method to quickly check the existence of LCO
close all
% CASE I: LCO with unstable BE
 A_1 = [-5 1 0;-9 0 1;-5 0 0];
 
 A_2 = [-5 1 0;12 0 1;-14 0 0];
 eig(A_1)
 eig(A_2)
 
 % 
 t_1 = [1.5:0.02:5];
 
 t_2 = [2:0.02:5];
 
 n_1 = length(t_1);
 n_2 = length(t_2);
 n_m = size(A_1,1);
 % 
 DET = zeros(n_1,n_2);
 
 for i =1:n_1
     for j=1:n_2
         temp = expm(A_1*t_1(i))*expm(A_2*t_2(j));
         DET(i,j) = det(temp - eye(n_m));
     end
 end
 
 X =repmat(t_1',1,n_2);
 Y = repmat(t_2,n_1,1);
 % 
 FIG1 = figure;
 plot3(X,Y,DET,'.');
 xlabel('t_1')
 ylabel('t_2')
 zlim([0 100])
%  any(DET>0)