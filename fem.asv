%% Setup
clear;
close all;
clear;
left = 0; % boundaries
right = 1;
m = 40; % number of points
M = zeros(m); % initialization
L = zeros(m);
K = zeros(m);
u0 = zeros([m,1]);
x = zeros([m,1]);
h = (right-left)/(m);


u_0 = 1; % amplitude?
k = 2*pi*2; % wave speed?
analycic = @(x,t) u_0*exp(1i*k*(x-t)); % might be the analythical solution


for i = 1:m
    x(i) = h*(i-1);
    u0(i) = real(analycic(x(i),0));
end
u1 = u0;
u2 = [u0, u0];


%% Assemble mass and stiffness matrix
mm = [1, 4, 1]; % weights for mass and stiffness matrix 
ll = [-1, 0, 1];
kk = [-1,2,-1];

for i = 1:m
    for j = 1:3
        M(i,mod(i+j-3,m)+1) = mm(j);
    end
    for j = 1:3
        L(i,mod(i+j-3,m)+1) = ll(j);
    end
    for j = 1:3
        K(i,mod(i+j-3,m)+1) = kk(j);
    end
end

M = M*(h/6);
L = L/2;
K = K/h;



%% Time step operator construncton
T = 0;
dt = 0.01;
%a = h/8;
a = 0;
%a = 0.000000000000000001*4.09;
%a = h/2000;
O = (M+dt/2*(L+a*K))\(M-dt/2*(L+a*K));% Crank nicolson
RK = -M\(L+a*K)*dt;
O2 = -M\(L+a*K)*2*dt; % leap frog
% Plot eigenvalued for operator
figure
ei = eig(RK);
plot(real(ei),imag(ei),'*');
maxdt = 2.7/max(imag(ei));% maxdt = 2.58/max(abs(ei));
disp(['Maximum RK4 timestep? - ' num2str(maxdt)])
pause
%% Time stepping
figure;
while T < 10
    T = T + dt;
   %u1 = u1 - dt*Central2*u1;
   %u1 = O*u1;
   
   %u1 = u1 + RK*u1;
   
   
   f1 =  RK * u1;
   f2 =  RK * (u1 + f1/2);
   f3 =  RK * (u1 + f2/2);
   f4 =  RK * (u1 + f3);
   u1 = u1 + (f1 + 2*f2 + 2*f3 + f4)/6;
   
   
   %{
   temp = u2(:,2);
   u2(:,2) = u2(:,1) + O2*u2(:,2); % leap frog
   u2(:,1) = temp;
   %}
   plot(x,real(analycic(x,T)),x,u1);
   title(T)
   pause(0.01)
end
figure;
plot(x,real(analycic(x,1)),x,u1);
legend('analytic','Numerical',"Location","best");
title(T)