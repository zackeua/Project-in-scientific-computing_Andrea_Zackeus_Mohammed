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
analycic = @(x,t) u_0*exp(1i*k*(x-t)); % might be the analythical solution

for i = 1:m
    x(i) = h*(i-1);
    u0(i) = analycic(x(i),0);
end


u1 = u0;

%% create differential operators
% weights for second and 4th order central difference
for degree = 2*(1:8)
    Central = zeros(m);
    w = weights(degree,1);
    for i = 1:m
        for j = 1:length(w)
            Central(i,mod(i+j-2-floor(length(w)/2),m)+1) = w(j);
        end
    end
    ei = eig(Central);
    figure
    plot(ei,'*');
    title(['Eigenvalues for h^{', num2str(degree),'}'])
    xlabel('Re(\lambda)')
    ylabel('Im(\lambda)')
end
