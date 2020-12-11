close all;

load('uniformStability');
C_uni = C_eff;

load('gaussLobattoStability');
C_gaussLobatto = C_eff;

load('gaussLobattoExactStability');
C_gaussLobattoExact = C_eff;


figure
plot(C_uni,'-*','LineWidth',10);
hold on;
plot(C_gaussLobatto,'-*','LineWidth',10);
plot(C_gaussLobattoExact,'-*','LineWidth',10);
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