% using 2D searching method to quickly check the existence of LCO
 close all
% CASE I: LCO with unstable BE
%  A_1 = [0 -1; 1 0];
 %
 A_1 = [1j 0; 0 1j];
%  A_2 = [-5 1 0;12 0 1;-14 0 0];
 eig(A_1)
%  eig(A_2)
 
 % 
 t_1 = [0:0.02:10];
 

 
 n_1 = length(t_1);
 
 n_m = size(A_1,1);
 % 
 DET = zeros(n_1,1);
 
 for i =1:n_1
  
         temp = expm(A_1*t_1(i));
         DET(i) = det(temp - eye(n_m)); 
 end
 
%  X =repmat(t_1',1,n_2);
%  Y = repmat(t_2,n_1,1);
 % 
 FIG1 = figure;
 plot(t_1,DET,'.');
 xlabel('t_1')
 ylabel('t_2')
 
 
%  zlim([0 100])
%  any(DET>0)

[v,d] = eig(expm(A_1*2*pi))
[v,d] = eig(expm(A_1*0.5*pi))
% CONCLUSION:

% This gives the imaginary eigenvalue's period with the det(K) fram.