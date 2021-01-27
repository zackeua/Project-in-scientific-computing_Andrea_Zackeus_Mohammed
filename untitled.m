close all;

%a = '1';
%a = '2';
%a = 10;
%a = 0.001;
degrees = 1:10;
bounds = [0,1];
n = 10;

scale = fem(0,1,1,bounds,n);

C_0 = fem(0,2,degrees,bounds,n);
C_0_001 = fem(0.001,2,degrees,bounds,n);
C_0_01 = fem(0.01,2,degrees,bounds,n);
C_h = fem('2',2,degrees,bounds,n);

C_h = C_h./scale;
C_0_01 = C_0_01./scale;
C_0_001 = C_0_001./scale;
C_0 = C_0./scale;

width = 5;

figure('Renderer', 'painters', 'Position', [10 10 1300 700])
plot(1+degrees,C_h,'-*','LineWidth',width);
hold on;
plot(1+degrees,C_0_01,'-*','LineWidth',width);
plot(1+degrees,C_0_001,'-*','LineWidth',width);
plot(1+degrees,C_0,'-*','LineWidth',width);
hold off;
xlim([1+degrees(1), 1+degrees(end)]);
legend('\alpha \propto (\Delta x_{min})^{2}', '\alpha = 0.01', '\alpha = 0.001', '\alpha = 0', 'Location', 'northeast');
xlabel('Order of accuracy')
ylabel('C_{eff}')
title('Rescaled C_{eff} for stabilized FEM with GL points and quadrature')
%title(['Rescaled C_{eff} for different numerical methods with \alpha = {', num2str(a),'}']);
%title(['Rescaled C_{eff} for different numerical methods with \alpha \propto \Delta x_{min}']);
set(gca,'FontSize',30)

figure('Renderer', 'painters', 'Position', [10 10 1300 700])
plot(1+degrees,1./C_h,'-*','LineWidth',width);
hold on;
plot(1+degrees,1./C_0_01,'-*','LineWidth',width);
plot(1+degrees,1./C_0_001,'-*','LineWidth',width);
plot(1+degrees,1./C_0,'-*','LineWidth',width);
hold off;
ylim([0,30])
xlim([1+degrees(1), 1+degrees(end)]);
legend('\alpha \propto (\Delta x_{min})^{2}', '\alpha = 0.01', '\alpha = 0.001', '\alpha = 0', 'Location', 'northeast');
xlabel('Order of accuracy');
ylabel('Reciprocal of C_{eff}')
title('Reciprocal rescaled C_{eff} for stabilized FEM with GL points and quadrature')
%title(['Reciprocal rescaled C_{eff} for different numerical methods with \alpha = {', num2str(a),'}'])
%title(['Reciprocal rescaled C_{eff} for different numerical methods with \alpha \propto \Delta x_{min}']);
set(gca,'FontSize',30)