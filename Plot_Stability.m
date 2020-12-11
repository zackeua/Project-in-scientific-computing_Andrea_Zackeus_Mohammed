% Computation of the A-Stability Region
% Borrowed from 
% http://www.sci.sdsu.edu/Faculty/Don.Short/math542/matlabcode.htm#A-STABILITY%20REGION%20CALCULATION
clear all;close all;

x=linspace(-3,1.5,200);
y=linspace(-3.5,3.5,200);

[X,Y]=meshgrid(x,y);
Z=X +Y*1i;
%Euler's Method
M=abs(1+Z);
contour(X,Y,M,[1,1],'r');
hold on
%Heun's Method
M=abs(1+Z+Z.^2/2);
contour(X,Y,M,[1,1],'g');
%RK4
M=abs(1+Z+Z.^2/2+Z.^3/6+Z.^4/24);
contour(X,Y,M,[1,1],'b');
grid on
axis equal
axis([-5 1 -3 3])
title('Runge-Kutta A-Stability Regions')
legend('Eulers Method', 'Heuns Method','RK4', 'Location', 'northwest')
xlabel('Re')
ylabel('Im')
set(gca,'FontSize',15);