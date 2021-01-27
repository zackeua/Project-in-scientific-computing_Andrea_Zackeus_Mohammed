close all;

%a = '1';
%a = '2';
%a = 10;
%a = 0.001;
degrees = 1:10;
bounds = [0,1];
n = 10;
C_uni = fem(a,1,degrees,bounds,n);

C_gaussLobatto = fem(a,2,degrees,bounds,n);

C_gaussLobattoExact = fem(a,3,degrees,bounds,n);

C_gaussLobatto = C_gaussLobatto./C_uni(1);
C_gaussLobattoExact = C_gaussLobattoExact./C_uni(1);
C_uni = C_uni./C_uni(1);

width = 5;

figure('Renderer', 'painters', 'Position', [10 10 1300 700])
plot(1+degrees,C_uni,'-*','LineWidth',width);
hold on;
plot(1+degrees,C_gaussLobatto,'-*','LineWidth',width);
plot(1+degrees,C_gaussLobattoExact,'-*','LineWidth',width);
hold off;
xlim([1+degrees(1), 1+degrees(end)]);
legend('FEM - evenly spaced points with exact integration', 'FEM - GL points with GL quadrature', 'FEM - GL points with exact integration', 'Location', 'northeast');
xlabel('Order of accuracy')
ylabel('C_{eff}')
title(['Rescaled C_{eff} for different numerical methods with \alpha = {', num2str(a),'}']);
title(['Rescaled C_{eff} for different numerical methods with \alpha \propto \Delta x_{min}']);
set(gca,'FontSize',30)

figure('Renderer', 'painters', 'Position', [10 10 1300 700])
plot(1+degrees, 1./C_uni,'-*','LineWidth',width);
hold on;
plot(1+degrees, 1./C_gaussLobatto,'-*','LineWidth',width);
plot(1+degrees, 1./C_gaussLobattoExact,'-*','LineWidth',width);
hold off;
xlim([1+degrees(1), 1+degrees(end)]);
legend('FEM - evenly spaced points with exact integration', 'FEM - GL points with GL quadrature', 'FEM - GL points with exact integration', 'Location', 'northwest');
xlabel('Order of accuracy');
ylabel('Reciprocal of C_{eff}')
%title(['Reciprocal rescaled C_{eff} for different numerical methods with \alpha = {', num2str(a),'}'])
title(['Reciprocal rescaled C_{eff} for different numerical methods with \alpha \propto \Delta x_{min}']);
set(gca,'FontSize',30)