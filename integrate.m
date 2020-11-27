function [M,L,K] = integrate(degree,h,n)
%  ha x vector ist?llet f?r h, h(i)=x(i)-x(i-1)

% Integrates "exactly" P(degree) elements with a uniform spatial step length h and n nodal
% points and assembles Mass (M), Load (L) and Stiffness (K) matrices.
% assembly of FEM matrixes for equidistant grid

%degree = 2; % degree of FEM basis functions
%h = 1; % space step length
phi = coeff(degree,h); % generate coefficients of basis functions
phiPrim = zeros(degree+1,degree);% vi forlorar en kolumn pga derivering
for i=1:degree+1 % differentiate basis functions
    phiPrim(i,:) = polyder(phi(i,:));
end

%n = degree+1; % for the size of the matrixes
m = zeros(degree+1);
l = zeros(degree+1);
%M = zeros((degree+1)+n*degree-1);% samma sak som 
M = zeros(degree*(n+1));% n+1:interval
L = zeros((degree+1)+n*degree-1);
K = zeros((degree+1)+n*degree-1);% lokal to global % slutade fraga har
a = 0;
b = degree*h;
for j = 1:degree+1
    for i = 1:degree+1
        m(i,j) = diff(polyval(polyint(conv(phi(j,:),phi(i,:))),[a,b]));
        l(i,j) = diff(polyval(polyint(conv(phiPrim(j,:),phi(i,:))),[a,b]));
        k(i,j) = diff(polyval(polyint(conv(phiPrim(j,:),phiPrim(i,:))),[a,b]));
    end
end

for i = 1:degree:length(M)-degree % Assemble FEM matrixes
    M(i:degree+i,i:degree+i) = M(i:degree+i,i:degree+i) + m;
    L(i:degree+i,i:degree+i) = L(i:degree+i,i:degree+i) + l;
    K(i:degree+i,i:degree+i) = K(i:degree+i,i:degree+i) + k;
end


% to fix last addition of local 
M(1,1) = M(1,1) + m(end,end);
L(1,1) = L(1,1) + l(end,end);
K(1,1) = K(1,1) + k(end,end);

M(end-degree+1:end,end-degree+1:end) = M(end-degree+1:end,end-degree+1:end) + m(1:degree,1:degree);
L(end-degree+1:end,end-degree+1:end) = L(end-degree+1:end,end-degree+1:end) + l(1:degree,1:degree);
K(end-degree+1:end,end-degree+1:end) = k(end-degree+1:end,end-degree+1:end) + k(1:degree,1:degree);

M(1,end-degree+1:end) = M(1,end-degree+1:end) + m(end,1:degree);
L(1,end-degree+1:end) = L(1,end-degree+1:end) + l(end,1:degree);
K(1,end-degree+1:end) = K(1,end-degree+1:end) + k(end,1:degree);

M(end-degree+1:end,1) = M(end-degree+1:end,1) + m(1:degree,end);
L(end-degree+1:end,1) = L(end-degree+1:end,1) + l(1:degree,end);
K(end-degree+1:end,1) = K(end-degree+1:end,1) + k(1:degree,end);



% print matrixes
%M % mass matrix
%L % phiPrim * phi matrix
%D % phiPRim * phiPrim matrix (stabilization)
end
