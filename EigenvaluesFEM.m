%% Setup
clear;
close all;
clc;
left = 0; % boundaries
right = 1;

%m = 60; % number of spatial points has to be dividable evenly by degree
n=10;

plotting = 1; % s??tt till 1 om du vill plotta

% Gauss-Lobatto points: v??lj m+1 punkter och ta bort sista
% [x,w]= legendre_gauss_lobatto(m+1);
% x= (right-left)/2 * x +(right +left)/2;
% w= w*((right-left)/2); % according to what Gunilla said
% x=flip(x);
% x = x(1:end-1);

%% analytic solution
u_0 = 1; % amplitude
k = 2*pi; % wave frequency
analytic = @(x,t) real(u_0*exp(1i*k*(x-t)));


%% Stability region of RK4 och RK1/explicit euler
xi = -3:0.01:0.5;
yi = -3:0.01:3;
[xi,yi] = meshgrid(xi,yi);
z = xi + 1i*yi;
stabRK4 = abs(1+z+1/2*z.^2+1/6*z.^3 + 1/24*z.^4);
stabRK1 = abs(1+z);

%m = (degree+1)+n*degree-1
for degree = 1:18
    % antal intervalinterval
    m=degree*(n+1);
    
    u0 = zeros([m,1]);
    x = zeros([m,1]);
    h = (right-left)/(m);
    
    %% Evenly spaced interpolation points or Gauss-Lobatto interpolation points
    for i = 1:m % evenly spaced points
        x(i) = h*(i-1);
    end
    
    for i = 1:m
        u0(i) = analytic(x(i),0);
    end
    
    %% Assemble mass and stiffness matrix
    % For evenly spaced points (using or not using Masslumping)
    [M,L,K] = integrate(degree,h,n);
    % For Gauss-Lobatto points and Gauss-Lobatto quadrature
% [M,L,K] = integrate2_GaussLobatto(degree,x,w);
    a = 0;
    %a = h/2000;
    RK = -M\(L+a*K);
    RK_Masslumping = -(eye(m).*sum(M))\(L+a*K);

    ei = eig(RK);
    ei_Masslumping = eig(RK_Masslumping);
    
    eimax = max(abs(ei));
    dtmax = 2.83/eimax;
    
    eimax = max(abs(ei_Masslumping));
    dtmax_Masslumping = 2.83/eimax;
    
    
    disp(['Order ', num2str(degree), ' biggest possible timestep ', num2str(dtmax)])
    disp(['Order ', num2str(degree), ' with masslumping biggest possible timestep ', num2str(dtmax_Masslumping)])
    
    
    if plotting == 1
        figure;
        contour(xi,yi,stabRK4,[1,1],'r')
        hold on;
        contour(xi,yi,stabRK1,[1,1])
        plot(ei*dtmax,'*b');
        plot(ei_Masslumping*dtmax_Masslumping,'+g');
        title(['Eigenvalues for P', num2str(degree),' elements'])
        xlabel('Re(\lambda)')
        ylabel('Im(\lambda)')
        ylim([-3 3]);
        xlim([-3 3]);
        axis equal;
        legend('RK4','explicit euler', '\lambda\cdot dt', '\lambda_{Masslumping}\cdot dt', 'Location','best')
        hold off;
    end
        
    %% w/o masslumping timestepping
    if plotting == 1
        u1 = u0;
        dt = dtmax*0.9;
        T = 0;
        while T < 30
           g1 = dt * RK * u1; %% minus ??r inlaggt i RK = -M\(L+a*K)
           g2 = dt * RK * (u1 + g1/2);
           g3 = dt * RK * (u1 + g2/2);
           g4 = dt * RK * (u1 + g3);
           u1 = u1 + (g1 + 2*g2 + 2*g3 + g4)/6;
           T = T + dt;
        end
    
        figure;
        plot([x; 1],[analytic(x,T); analytic(x(1),T)],[x; 1], [u1; u1(1)]);
        legend('analytic','RK4',"Location","best");
        title(['T = ', num2str(T), ', with P', num2str(degree), ' elements and timestep: ', num2str(dt)])
    end
    
    %% MAss lumping time stepping
    if plotting == 1
        u1 = u0;
        dt =  dtmax_Masslumping*0.9;
        T = 0;
        while T < 30
           g1 = dt * RK_Masslumping * u1; %% minus ??r inlaggt i RK = -M\(L+a*K)
           g2 = dt * RK_Masslumping * (u1 + g1/2);
           g3 = dt * RK_Masslumping * (u1 + g2/2);
           g4 = dt * RK_Masslumping * (u1 + g3);
           u1 = u1 + (g1 + 2*g2 + 2*g3 + g4)/6;
           T = T + dt;
        end
    
        figure;
        plot([x; 1],[analytic(x,T); analytic(x(1),T)],[x; 1], [u1; u1(1)]);
        legend('analytic','RK4',"Location","best");
        title(['T = ', num2str(T), ', with P', num2str(degree), ' elements and timestep: ', num2str(dt), ' with masslumping'])
    end
    
    
    
    %e = analytic(x,T) - u1;
    %E = norm(e'*e);
  
end
%% The rescaled efficiency number
degree = 1:6;
C_eff= zeros(1,6);
for i= 1:6
    % the wave number
    c= 1; % not sure if it should be equal to one??
%     c=2*pi;
    % the width of the interval
    width= 1;
    % total number of degrees of freedom % m/degree-1=antalet elements
    % N_dof = ( m/i-1)*i;
    N_dof = m*i;
    C_eff(i) = (sqrt(3)*c*dt*N_dof)/width;
end
figure;
plot(degree,C_eff,'-*','MarkerIndices',1:length(C_eff) );
xlabel('Polynomial degree');
ylabel('Rescaled C_{eff}');
title('Rescaled efficiency number as a function of the polynomial degree  ')