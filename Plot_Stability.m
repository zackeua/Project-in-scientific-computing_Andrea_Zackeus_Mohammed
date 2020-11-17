% Computation of the A-Stability Region
% Borrowed from 
% http://www.sci.sdsu.edu/Faculty/Don.Short/math542/matlabcode.htm#A-STABILITY%20REGION%20CALCULATION
clear all;close all;

x=linspace(-3,1.5,200);
y=linspace(-3.5,3.5,200);

[X,Y]=meshgrid(x,y);
Z=X +Y*i;
%Euler's Method
M=abs(1+Z);
[c,h]=contour(X,Y,M,[1,1]);
set(h,'linewidth',2,'edgecolor','b')
hold on
%Heun's Method
M=abs(1+Z+Z.^2/2);
[c,h]=contour(X,Y,M,[1,1]);
set(h,'linewidth',2,'edgecolor','g')
%RK4
M=abs(1+Z+Z.^2/2+Z.^3/6+Z.^4/24);
[c,h]=contour(X,Y,M,[1,1]);
set(h,'linewidth',2,'edgecolor','r','lineStyle',':')

grid on
axis equal
axis([-5 1 -3 3])
title('Runge-Kutta A-Stability Regions')
legend('Eulers Method', 'Heuns Method','RK4')
xlabel('Re')
ylabel('Im')