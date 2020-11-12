%% Setup
clear;
close all;
clear;
left = 0; % boundaries
right = 1;
m = 20; % number of points
Central2 = zeros(m); % initialization
Central4 = zeros(m);
u0 = zeros([m,1]);
x = zeros([m,1]);
h = (right-left)/(m);


u_0 = 1; % amplitude?
k = 2*pi; % wave speed?
analycic = @(x,t) u_0*exp(1i*k*(x-t)); % might be the analythical solution


for i = 1:m
    x(i) = h*(i-1);
    u0(i) = analycic(x(i),0);
end


u1 = u0;
u2 = u0;
u3 = u0;
u4 = u0;
u5 = [u0,u0];
u6 = u5;
%plot(x,real(u0))

%% create differential operators
% weights for second and 4th order central difference
T1 = [-1, 0, 1];
T2 = [1, -8, 0, 8, -1];

for i = 1:m
    for j = 1:3
        Central2(i,mod(i+j-3,m)+1) = T1(j);
    end
    for j = 1:5
        Central4(i,mod(i+j-4,m)+1) = T2(j);
    end 
end
Central2
Central4
Central2 = Central2/(2*h);
Central4 = Central4/(12*h);

%% eigenvalue plotting

ei2 = eig(Central2*h);
ei4 = eig(Central4*h);

figure
plot(real(ei2),imag(ei2),'+',real(ei4),imag(ei4),'*');
title('Eigenvalues')
xlabel('Re(\lambda)')
ylabel('Im(\lambda)')
legend('h^2','h^4','location','best')

%% Time stepping
% euler forward
% (u_(t+1) - u_t)/dt +O*u_t = 0 
% (u_(t+1) - u_t)/dt = -O*u_t
% u_(t+1) - u_t = -dt*O*u_t
%  u_(t+1) = u_t -dt*O*u_t

% du/dt = F(t,u) = -O*u

% f1 = F(t,u)
% f2 = F(t,u+f1/2)
%...

T = 0;
dt = 0.01;

while T < 1
   
   f1 = -dt * Central2 * u1;
   f2 = -dt * Central2 * (u1 + f1/2);
   f3 = -dt * Central2 * (u1 + f2/2);
   f4 = -dt * Central2 * (u1 + f3);
   u1 = u1 + (f1 + 2*f2 + 2*f3 + f4)/6;
   
   u2 = u2 - dt*Central2*u2;
   
   g1 = dt * Central4 * u3;
   g2 = dt * Central4 * (u3 + g1/2);
   g3 = dt * Central4 * (u3 + g2/2);
   g4 = dt * Central4 * (u3 + g3);
   u3 = u3 + (g1 + 2*g2 + 2*g3 + g4)/6;
   
   u4 = u4 - dt*Central4*u4;
   temp = u5(:,2);
   u5(:,2) = u5(:,1) - 2*dt*Central2*u5(:,2);
   u5(:,1) = temp;
   
   temp = u6(:,2);
   u6(:,2) = u6(:,1) - 2*dt*Central4*u6(:,2);
   u6(:,1) = temp;
   
   T = T + dt;
end
figure;
plot(x,analycic(x,1),x,u1,x,u2,x,u3,x,u4,x,u5(:,2),x,u6(:,2));
legend('analytic','RK4 - h^2','Euler - h^2','RK4 - h^4','Euler - h^4','leapfrog - h^2','leapfrog - h^4',"Location","best");
title(T)