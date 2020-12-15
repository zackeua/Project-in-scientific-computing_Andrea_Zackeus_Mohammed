function [C_eff,L] = fem(a,choice,degrees,bounds,intervals)
% a is stabilisation term in the equation u_t + u_x = a*u_xx
% choice represents, 1=uniform, 2=gauss lobatto quadrature, 3=GL with exact
% integration
% degrees is a vector of degrees of intrest
% bounds is [left, right] bounds of domain, i.e. [0,1]
% intervals is the number of intervals to split the domain into
%setup
left = bounds(1); % boundaries
right = bounds(2);

%m = 60; % number of spatial points has to be dividable evenly by degree
n=intervals;

C_eff = [];

plotting = 0; % s??tt till 1 om du vill plotta

plot_eigenvalues = 1; % välj vad du vill plotta och skriva ut
plot_C_eff = 0;
disp_max_timesteps = 0;

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
for degree = degrees
    % antal intervalinterval
    m=degree*(n+1);
    
    %% Assemble mass and stiffness matrix
    % For evenly spaced points (using or not using Masslumping)
    [M,L,K,X] = MatrixAssembler(degree,n,choice,bounds);
    u0 = analytic(X,0);
    % For Gauss-Lobatto points and Gauss-Lobatto quadrature
    % [M,L,K] = integrate2_GaussLobatto(degree,x,w);
    %a = 0;
    %h = X(2)-X(1);
    %L(1:degree,1:degree)
    %K(1:degree,1:degree)
    h_vec = [X(2:end); right] - X;
    %h = min(h_vec);
    
    %a = h; %h*h;
    RK = M\(-a*K-L);

    ei = eig(RK);
     
    eimax = max(abs(ei));
    dtmax = 2.5/eimax;
    
    C_eff = [C_eff sqrt(3)*dtmax*m]; % calculate next C_eff number only
    
    if disp_max_timesteps == 1
        disp(['Order ', num2str(degree), ' biggest possible timestep ', num2str(dtmax)])
    end
    
    if plot_eigenvalues == 1
        figure;
        contour(xi,yi,stabRK4,[1,1],'r')
        hold on;
        contour(xi,yi,stabRK1,[1,1])
        plot(ei*dtmax,'*b');
        title(['Eigenvalues for P', num2str(degree),' elements'])
        xlabel('Re(\lambda)')
        ylabel('Im(\lambda)')
        ylim([-3 3]);
        xlim([-3 3]);
        %axis equal;
        legend('RK4','explicit euler', '\lambda\cdot dt', 'Location','best')
        hold off;
    end
    %% w/o masslumping timestepping
    if plotting == 1
        %pause;
        u1 = u0;
        dt = dtmax*0.9*0.1;
        T = 0;
        figure;
        while T < 1
           g1 = dt * RK * u1; %% minus ??r inlaggt i RK = -M\(L+a*K)
           g2 = dt * RK * (u1 + g1/2);
           g3 = dt * RK * (u1 + g2/2);
           g4 = dt * RK * (u1 + g3);
           u1 = u1 + (g1 + 2*g2 + 2*g3 + g4)/6;
           T = T + dt;
           plot([X; right],[analytic(X,T); analytic(X(1),T)],[X; right], [u1; u1(1)]);
        legend('analytic','RK4',"Location","best");
        title(['T = ', num2str(T), ', with P', num2str(degree), ' elements and timestep: ', num2str(dt)])
        pause(0.01)
        end
        
        plot([X; right],[analytic(X,T); analytic(X(1),T)],[X; right], [u1; u1(1)]);
        legend('analytic','RK4',"Location","best");
        title(['T = ', num2str(T), ', with P', num2str(degree), ' elements and timestep: ', num2str(dt)])
    end
    
end
%% The rescaled efficiency number
if plot_C_eff == 1
    figure;
    plot(degrees,C_eff,'-*','MarkerIndices',1:length(C_eff));
    xlabel('Polynomial degree');
    ylabel('Rescaled C_{eff}');
    title('Rescaled efficiency number as a function of the polynomial degree');

    figure;
    plot(degrees,1./C_eff,'-*','MarkerIndices',1:length(C_eff));
    xlabel('Polynomial degree');
    ylabel('1/Rescaled C_{eff}');
    title('Rescaled efficiency number as a function of the polynomial degree');
end

end