Lie = @(x,y) x*y -y*x;
A1 = [-3.2 -1 0; 25.61 0 -1; -75.03 0 0];
A2 = [-1 -1 0; 1.28 0 -1; -0.624 0 0];
Lie(A1,A2)
%
A1 = [-1 -1;1 -1]; A2 = [-1 -2; 0.5 -1];
L1 =Lie(A1,A2)
L21 = Lie(A1,L1) 
L22= Lie(A2,L1)
L2 = L21 + L22
L3 =Lie(A2,L21)

