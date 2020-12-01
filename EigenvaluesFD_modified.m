%% Setup
% compy of EigenvaluesFD
clear;
close all;
clc;
left = 0; % boundaries
right = 1;
n=10;

plot_CFL=1;
CFL=[];

%m = 60; % number of points has to be at least degree*2+1
% m=degree*(n+1);
% Central = zeros(m); % initialization
% u0 = zeros([m,1]);
% x = zeros([m,1]);
% h = (right-left)/(m);


u_0 = 1; % amplitude
k = 2*pi; % wave frequency
analytic = @(x,t) real(u_0*exp(1i*k*(x-t)));






%%
%f1 = @ (x) ();
%function semicirc = semicircle(x, y)
r = 2.58;
phi = linspace(-pi,0,1000);
c2 = r*cos(phi); c2 = [c2,c2(1)];
c1 = r*sin(phi); c1 = [c1,c1(1)];
%patch(x,y,'y');
%axis equal;
xi = -3:0.01:0.5;
yi = -3:0.01:3;
[xi,yi] = meshgrid(xi,yi);
z = xi + 1i*yi;
stabRK4 = abs(1+z+1/2*z.^2+1/6*z.^3 + 1/24*z.^4);
stabRK1 = abs(1+z);


%% create differential operators
%weights for second and 4th order central difference
for degree = 2*(1:8)
    m=degree*(n+1);
    Central = zeros(m); % initialization
    h = (right-left)/(m);
    x = zeros([m,1]);
    u0 = zeros([m,1]);
    for i = 1:m
        x(i) = h*(i-1);
        u0(i) = analytic(x(i),0);
    end
    w = -weights(degree,1)/h;
    for i = 1:m
        for j = 1:length(w)
            Central(i,mod(i+j-2-floor(length(w)/2),m)+1) = w(j);
        end
    end
    ei = eig(Central);
    eimax = max(abs(ei));
    dtmax = 2.83/eimax;
    disp(['Order ', num2str(degree), ' biggest possible timestep ', num2str(dtmax)])
    
    CFL=[CFL dtmax/h];
    
    
    figure
    %patch(c1,c2,'y');
    contour(xi,yi,stabRK4,[1,1],'r')
    hold on;
    contour(xi,yi,stabRK1,[1,1])
    plot(ei*dtmax,'*b');
    title(['Eigenvalues for finite differece of order ', num2str(degree)])
    xlabel('Re(\lambda)')
    ylabel('Im(\lambda)')
    legend('RK4','explicit euler', '\lambda\cdot dt','Location','best')
    hold off;
    
    
    u1 = u0;
    dt = dtmax*0.9;
    T = 0;
    while T < 30
       g1 = dt * Central * u1;%% la minus i Central = -Central
       g2 = dt * Central * (u1 + g1/2);
       g3 = dt * Central * (u1 + g2/2);
       g4 = dt * Central * (u1 + g3);
       u1 = u1 + (g1 + 2*g2 + 2*g3 + g4)/6;
       T = T + dt;
    end
    figure;
    plot([x; 1],[analytic(x,T); analytic(x(1),T)],[x; 1], [u1; u1(1)]);
    legend('analytic','RK4',"Location","best");
    title(['T = ', num2str(T), ', with order: ', num2str(degree), ' and timestep: ', num2str(dt)])
end
degree=  2*(1:8);
 %% The CFL number
    if plot_CFL == 1
        figure;
        plot(degree,CFL,'-*','MarkerIndices',1:length(CFL));
        xlabel('Polynomial degree');
        ylabel('CFL number');
        title('CFL number as a function of the polynomial degree');
    end
