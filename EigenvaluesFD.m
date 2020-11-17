%% Setup
clear;
close all;
clear;
left = 0; % boundaries
right = 1;
m = 40; % number of points
Central = zeros(m); % initialization
u0 = zeros([m,1]);
x = zeros([m,1]);
h = (right-left)/(m);


u_0 = 1; % amplitude?
k = 2*pi; % wave speed?
analycic = @(x,t) real(u_0*exp(1i*k*(x-t))); % might be the analythical solution

for i = 1:m
    x(i) = h*(i-1);
    u0(i) = analycic(x(i),0);
end




%%
%f1 = @ (x) ();
%function semicirc = semicircle(x, y)
r = 2.58;
phi = linspace(-pi,0,1000);
c2 = r*cos(phi); c2 = [c2,c2(1)];
c1 = r*sin(phi); c1 = [c1,c1(1)];
%patch(x,y,'y');
%axis equal;

%% create differential operators
%weights for second and 4th order central difference
for degree = 2*(1:8)
    Central = zeros(m);
    w = weights(degree,1);
    for i = 1:m
        for j = 1:length(w)
            Central(i,mod(i+j-2-floor(length(w)/2),m)+1) = w(j);
        end
    end
    ei = eig(-Central);
    eimax = max(abs(ei));
    dtmax = r/eimax;
    maxReEi = max(real(ei));
    if (maxReEi*dtmax >= eps(1.))
        if (maxReEi ~= 0)
            dtmax = eps(1.)/maxReEi;
        end
    end
    disp(['Order ', num2str(degree), ' biggest possible timestep ', num2str(dtmax)])
%     
%     figure
%     patch(c1,c2,'y');
%     hold on;
%     plot(ei,'*');
%     title(['Eigenvalues for h^{', num2str(degree),'}'])
%     xlabel('Re(\lambda)')
%     ylabel('Im(\lambda)')
%     hold off;
%     
    
    u1 = u0;
    dt = dtmax*0.9/2;
    T = 0;
    while T < 1
       g1 = dt * Central * u1;
       g2 = dt * Central * (u1 + g1/2);
       g3 = dt * Central * (u1 + g2/2);
       g4 = dt * Central * (u1 + g3);
       u1 = u1 + (g1 + 2*g2 + 2*g3 + g4)/6;
       
       T = T + dt;
    end
    figure;
    %plot(x,analycic(x,1),x,u1,x,u2,x,u3,x,u4,x,u5(:,2),x,u6(:,2));
    %legend('analytic','RK4 - h^2','Euler - h^2','RK4 - h^4','Euler - h^4','leapfrog - h^2','leapfrog - h^4',"Location","best");
    plot(x,analycic(x,1),x,u1);
    legend('analytic','RK4',"Location","best");
    
    
    title(['T = ', num2str(T), ', with order: ', num2str(degree), ' and timestep: ', num2str(dt)])
end
