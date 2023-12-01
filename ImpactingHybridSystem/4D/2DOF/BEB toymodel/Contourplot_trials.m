% check the contour plot of the 3D
clc
clear 
close all
[X,Y,Z] = meshgrid(-2:.2:2);
V = X.*exp(-X.^2-Y.^2-Z.^2);

xslice = [-1.2,0.8,2];   
yslice = [];
zslice = [];
figure(1)
plot3(1,1,1,'r*')
hold on
contourslice(X,Y,Z,V,xslice,yslice,zslice)
view(3)
grid on