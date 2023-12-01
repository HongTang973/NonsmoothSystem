close all
clear 
% check the examples in book
N1=[-1 1;-1 0];
E=[3;-4];
C=[1;0];
N2=N1+E*C';
[V1,D1]=eig(N1)
[V2,D2]=eig(N2)