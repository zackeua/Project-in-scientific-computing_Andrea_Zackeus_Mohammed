close all;

a = '1';
degrees = 1:10;
bounds = [0,1];
n = 10;
C_uni = fem2(a,1,degrees,bounds,n);

C_gaussLobatto = fem2(a,2,degrees,bounds,n);

C_gaussLobattoExact = fem2(a,3,degrees,bounds,n);


figure('Renderer', 'painters', 'Position', [10 10 900 600])
plot(degrees,C_uni,'-*','LineWidth',10);
hold on;
plot(degrees,C_gaussLobatto,'-*','LineWidth',10);
plot(degrees,C_gaussLobattoExact,'-*','LineWidth',10);
hold off;
xlim([degrees(1), degrees(end)]);
legend('FEM - evenly spaced points with exact integration', 'FEM - GL points with GL quadrature', 'FEM - GL points with exact integration', 'Location', 'northeast');
xlabel('Degree')
ylabel('C_{eff}')
title('Rescaled C_{eff} for different numerical methods')
%title({'Rescaled efficiency number as a function of the','polynomial degree for FEM and', 'order of the FDM central scheme'});
set(gca,'FontSize',30)

figure('Renderer', 'painters', 'Position', [10 10 1100 600])
plot(1./C_uni,'-*','LineWidth',5);
hold on;
plot(1./C_gaussLobatto,'-*','LineWidth',5);
plot(1./C_gaussLobattoExact,'-*','LineWidth',5);
hold off;
legend('FEM - evenly spaced points with exact integration', 'FEM - GL points with GL quadrature', 'FEM - GL points with exact integration', 'Location', 'northwest');
xlabel('Degree')
ylabel('Reciprocal of C_{eff}')
title('Reciprocal rescaled C_{eff} for different numerical methods')
%title({'Reciprocal of rescaled efficiency number as a function of the','polynomial degree for FEM and', 'order of the FDM central scheme'});
set(gca,'FontSize',30)