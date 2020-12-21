close all;

a = '2';
degrees = 1:9;
bounds = [0,1];
n = 10;
C_uni = fem2(a,1,degrees,bounds,n);

C_gaussLobatto = fem2(a,2,degrees,bounds,n);

C_gaussLobattoExact = fem2(a,3,degrees,bounds,n);


figure
plot(degrees,C_uni,'-*','LineWidth',10);
hold on;
plot(degrees,C_gaussLobatto,'-*','LineWidth',10);
plot(degrees,C_gaussLobattoExact,'-*','LineWidth',10);
hold off;
legend('Lagrangian FEM - evenly spaced points with exact integration', 'Lagrangian FEM - GL points with GL quadrature', 'Lagrangian FEM - GL points with exact integration');
xlabel('Degree (FEM)')
ylabel('C_{eff}')
title('Rescaled efficiency number for different numerical methods')
%title({'Rescaled efficiency number as a function of the','polynomial degree for FEM and', 'order of the FDM central scheme'});
set(gca,'FontSize',30)

figure
plot(1./C_uni,'-*','LineWidth',10);
hold on;
plot(1./C_gaussLobatto,'-*','LineWidth',10);
plot(1./C_gaussLobattoExact,'-*','LineWidth',10);
hold off;
legend('Lagrangian FEM - evenly spaced points with exact integration', 'Lagrangian FEM - GL points with GL quadrature', 'Lagrangian FEM - GL points with exact integration', 'Location', 'northwest');
xlabel('Degree (FEM)')
ylabel('Reciprocal of C_{eff}')
title('Reciprocal rescaled efficiency number for different numerical methods')
%title({'Reciprocal of rescaled efficiency number as a function of the','polynomial degree for FEM and', 'order of the FDM central scheme'});
set(gca,'FontSize',30)