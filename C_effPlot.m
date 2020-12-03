close all;

load('uniform');
C_uni = C_eff;

load('gaussLobatto');
C_gaussLobatto = C_eff;

load('gaussLobattoExact');
C_gaussLobattoExact = C_eff;

load('fdm');


figure
plot(C_uni,'-*');
hold on;
plot(C_gaussLobatto,'-*');
plot(C_gaussLobattoExact,'-*');
plot(2*(1:4),CFL,'-*');
hold off;
legend('Uniform', 'Gauss Lobatto', 'Gauss Lobatto exact integration', 'FDM');
xlabel('Degree')
ylabel('C_{eff}')
title('Rescaled efficiency number as a function of the polynomial degree');


figure
plot(1./C_uni,'-*');
hold on;
plot(1./C_gaussLobatto,'-*');
plot(1./C_gaussLobattoExact,'-*');
plot(2*(1:4),1./CFL,'-*');
hold off;
legend('Uniform', 'Gauss Lobatto', 'Gauss Lobatto exact integration', 'FDM');
xlabel('Degree')
ylabel('1/C_{eff}')
title('Rescaled efficiency number as a function of the polynomial degree');
