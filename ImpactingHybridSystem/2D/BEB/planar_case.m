% planar system 16 kinds of orbits
% trace(A)： > < 0
% r: < > 1
% r*e^(pi/w) < >1
% node/focus

% % Detecting the existence of LCO when boundary equilibrium transferring to
% virtual Equilibrium for general n dimension impacting hybrid system 
% limitations: single impact; single discontinuity surface

% KEY INPUTs:
% A : the linear region system's matrix 
% R : the reset map matrix, R(x) = R*x;
% C : the selecting vector to define the discontinuity surface C*x=0

clc
clear
close all
% define the linear part
% As shown in the Exapmle 2.1, 
xi = -0.1; r = 0.5;
A  =  [0 1;-1 -2*xi];
MW = [0 0;0 -(1+r)];
R  = eye(size(A,1)) + MW;
C  = [1,0];