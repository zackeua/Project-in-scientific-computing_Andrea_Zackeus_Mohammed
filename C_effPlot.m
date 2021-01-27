close all;

a = 0;
degrees = 1:10;
bounds = [0,1];
n = 10;
C_uni = fem(a,1,degrees,bounds,n);

C_gaussLobatto = fem(a,2,degrees,bounds,n);


C_gaussLobattoExact = fem(a,3,degrees,bounds,n);


load('fdm');

C_dfem = [1.4200, 1.0800, 0.8845, 0.7684, 0.6766, 0.6120, 0.5592, 0.5142, 0.4788, 0.4501];


%order    1,      2,      3,      4,      5,      6,      7,      8,      9,      10
%dt       0.0710, 0.0360, 0.0221, 0.0154, 0.0113, 0.0087, 0.0070, 0.0057, 0.0048, 0.0041
%C_eff?   0.1420, 0.1080, 0.0884, 0.0768, 0.0677, 0.0612, 0.0559, 0.0514, 0.0479, 0.0450
%CFL      0.71,   0.72,   0.80,   0.89,   0.96,   1.03,   1.09,   1.14,   1.19,   1.24 
%         1.4200, 1.0800, 0.8845, 0.7684, 0.6766, 0.6120, 0.5592, 0.5142, 0.4788, 0.4501

C_gaussLobatto = C_gaussLobatto./C_uni(1);
C_gaussLobattoExact = C_gaussLobattoExact./C_uni(1);
CFL = CFL./C_uni(1);
C_dfem = C_dfem./C_uni(1);
C_uni = C_uni./C_uni(1);

width = 5;

figure('Renderer', 'painters', 'Position', [10 10 1300 700])
plot(1+degrees,C_uni,'-*','LineWidth',width);
hold on;
plot(1+degrees,C_gaussLobatto,'-*','LineWidth',width);
plot(1+degrees,C_gaussLobattoExact,'-*','LineWidth',width);
plot(1+degrees,C_dfem,'-*', 'LineWidth',width)
plot(2*(1:5),CFL,'-*','LineWidth',width);
hold off;
xlim([degrees(1)+1, degrees(end)+1]);
legend('FEM - evenly spaced points with exact integration', 'FEM - GL points with GL quadrature', 'FEM - GL points with exact integration', 'DGFEM - GL points with upwind flux for Jacobi polynomials', 'FDM - with central difference scheme', 'Location', 'northeast');
xlabel('Order of accuracy')
ylabel('C_{eff}')
title('Rescaled C_{eff} for different numerical methods')
%title({'Rescaled efficiency number as a function of the','polynomial degree for FEM and', 'order of the FDM central scheme'});
set(gca,'FontSize',30)

figure('Renderer', 'painters', 'Position', [10 10 1300 700])
plot(1+degrees, 1./C_uni,'-*','LineWidth',width);
hold on;
plot(1+degrees, 1./C_gaussLobatto,'-*','LineWidth',width);
plot(1+degrees, 1./C_gaussLobattoExact,'-*','LineWidth',width);
plot(1+degrees, 1./C_dfem,'-*', 'LineWidth',width)
plot(2*(1:5), 1./CFL,'-*','LineWidth',width);
hold off;
xlim([degrees(1)+1, degrees(end)+1]);
legend('FEM - evenly spaced points with exact integration', 'FEM - GL points with GL quadrature', 'FEM - GL points with exact integration', 'DGFEM - GL points with upwind flux for Jacobi polynomials', 'FDM - with central difference scheme', 'Location', 'northwest');
xlabel('Order of accuracy')
ylabel('Reciprocal of C_{eff}')
title('Reciprocal rescaled C_{eff} for different numerical methods')
%title({'Reciprocal of rescaled efficiency number as a function of the','polynomial degree for FEM and', 'order of the FDM central scheme'});
set(gca,'FontSize',30)

