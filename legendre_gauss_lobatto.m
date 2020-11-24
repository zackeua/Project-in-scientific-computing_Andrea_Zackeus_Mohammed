function [x,w] = legendre_gauss_lobatto (N)
% Nodes and weights for N-point Legendre Gauss-Lobatto quadrature
n = [0:N-3]'; beta = (n+1)./sqrt((2*n+1).*(2*n+3));
J_Nm1 = diag(beta,-1)+diag(beta,1);
gamma = (J_Nm1+eye(N-1))\[zeros(N-2,1); 1];
mu = (J_Nm1-eye(N-1))\[zeros(N-2,1); 1];
beta_Nm1 = sqrt(2/(gamma(N-1)-mu(N-1)));
alpha_N = 1+mu(N-1)*2/(gamma(N-1)-mu(N-1));
Jtilde_N = [J_Nm1 beta_Nm1*[zeros(N-2,1); 1]; beta_Nm1*...
[zeros(N-2,1); 1]' alpha_N]; % Symmetric tridiagonal matrix
[V,D] = eig(Jtilde_N); % Solve eigenvalue problem
x = diag(D); [x,i] = sort(x); % Gauss-Lobatto quadrature nodes
x = flipud(x); w = 2*V(1,i).^2'; % Gauss-Lobatto quadrature weights