clear;
close all;
clc;
ps = 2*(1:5);
CFL = [];
for degree = ps
    w = weights(degree,1);
    u = -degree/2:degree/2;

    deg = 0:0.01:2*pi;
    t = [];
    for angle = deg
        t = [t, abs(w*exp(1i*u*angle)')];
    end
    
    CFL = [CFL, sqrt(8)/max(t)];

    
    plot(deg, t);
    title('|t(\xi)|')
    xlabel('\xi');
    ylabel('|t(\xi)|')
    
    disp(['CFL = ' num2str(CFL(end)), ' for p = ', num2str(degree)]);
end

plot(ps, CFL, '*-')
save('fdm','CFL')