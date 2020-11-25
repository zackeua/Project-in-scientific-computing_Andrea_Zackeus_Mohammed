function [M,L,K] = integrate2_hermite(degree,x)
% Integrates "exactly" P(degree) elements with a uniform spatial step length h and n nodal
% points and assembles Mass (M), Load (L) and Stiffness (K) matrices.
% assembly of FEM matrixes for equidistant grid

%degree = 2; % degree of FEM basis functions
%h = 1; % space step length
pts = length(x);
%n = degree+1; % for the size of the matrixes
m = zeros(degree+1);
l = zeros(degree+1);
M = zeros(pts);
L = zeros(pts);
K = zeros(pts);% lokal to global % slutade fraga har


for i = 1:degree:length(M)-degree % Assemble FEM matrixes
    %mo
    %phi = coeff2(degree,x(i:degree+i)); % generate coefficients of basis functions
    phi = hermite_rec(degree);
    %
    phiPrim = zeros(degree+1,degree);% vi forlorar en kolumn pga derivering
%mo%     for j=1:degree+1 % differentiate basis functions
    for j=1:degree-1 % differentiate basis functions
        phiPrim(j,:) = polyder(phi(j,:));
    end
    a = x(i);
    b = x(i+degree);
   %mo% for j = 1:degree+1
   for j = 1:degree
      %mo%  for c = 1:degree+1
      c=1;
      while  c <= degree && (mod ( degree , 2 ) ~= 0)%mo%
          m(c,j) = diff(polyval(polyint(conv(phi(j,:),phi(c,:))),[a,b]));
          l(c,j) = diff(polyval(polyint(conv(phiPrim(j,:),phi(c,:))),[a,b]));
          k(c,j) = diff(polyval(polyint(conv(phiPrim(j,:),phiPrim(c,:))),[a,b]));
          c=c+1; %mo%
      end
   end
    k
    M(i:degree+i,i:degree+i) = M(i:degree+i,i:degree+i) + m;
    L(i:degree+i,i:degree+i) = L(i:degree+i,i:degree+i) + l;
    K(i:degree+i,i:degree+i) = K(i:degree+i,i:degree+i) + k;
end


% mo
%phi = coeff2(degree,[x(i+degree:end); 1]); % generate coefficients of basis functions
phi = hermite_rec(degree);
%
phiPrim = zeros(degree+1,degree);% vi forlorar en kolumn pga derivering
% mo% for j=1:degree+1 % differentiate basis functions
for j=1:degree
    phiPrim(j,:) = polyder(phi(j,:));
end
a = x(i+degree);
b = 1;
%mo% for j = 1:degree+1
%     for c = 1:degree+1
for j = 1:degree
    for c = 1:degree
        m(c,j) = diff(polyval(polyint(conv(phi(j,:),phi(c,:))),[a,b]));
        l(c,j) = diff(polyval(polyint(conv(phiPrim(j,:),phi(c,:))),[a,b]));
        k(c,j) = diff(polyval(polyint(conv(phiPrim(j,:),phiPrim(c,:))),[a,b]));
    end
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
