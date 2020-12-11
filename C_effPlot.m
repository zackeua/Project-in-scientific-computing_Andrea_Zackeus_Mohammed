close all;

load('uniform');
C_uni = C_eff;

load('gaussLobatto');
C_gaussLobatto = C_eff;

load('gaussLobattoExact');
C_gaussLobattoExact = C_eff;

load('fdm');


figure
plot(C_uni,'-*','LineWidth',10);
hold on;
plot(C_gaussLobatto,'-*','LineWidth',10);
plot(C_gaussLobattoExact,'-*','LineWidth',10);
plot(2*(1:4),CFL,'-*','LineWidth',10);
hold off;
legend('Lagrangian FEM - evenly spaced points with exact integration', 'Lagrangian FEM - GL points with GL quadrature', 'Lagrangian FEM - GL points with exact integration', 'FDM - with central difference scheme');
xlabel({'Degree (FEM)', 'Order (FDM)'})
ylabel('C_{eff}')
title('Rescaled efficiency number for different numerical methods')
%title({'Rescaled efficiency number as a function of the','polynomial degree for FEM and', 'order of the FDM central scheme'});
set(gca,'FontSize',30)

figure
plot(1./C_uni,'-*','LineWidth',10);
hold on;
plot(1./C_gaussLobatto,'-*','LineWidth',10);
plot(1./C_gaussLobattoExact,'-*','LineWidth',10);
plot(2*(1:4),1./CFL,'-*','LineWidth',10);
hold off;
legend('Lagrangian FEM - evenly spaced points with exact integration', 'Lagrangian FEM - GL points with GL quadrature', 'Lagrangian FEM - GL points with exact integration', 'FDM - with central difference scheme', 'Location', 'northwest');
xlabel({'Degree (FEM)','Order (FDM)'})
ylabel('Reciprocal of C_{eff}')
title('Reciprocal rescaled efficiency number for different numerical methods')
%title({'Reciprocal of rescaled efficiency number as a function of the','polynomial degree for FEM and', 'order of the FDM central scheme'});
set(gca,'FontSize',30)