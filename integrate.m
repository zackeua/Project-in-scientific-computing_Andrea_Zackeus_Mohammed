clear;
clc;

% assembly of FEM matrixes for equidistant grid

degree = 2; % degree of FEM basis functions
h = 1; % space step length
phi = coeff(degree,h); % generate coefficients of basis functions
phiPrim = zeros(degree+1,degree);
for i=1:degree+1 % differentiate basis functions
    phiPrim(i,:) = polyder(phi(i,:));
end

n = degree+1; % for the size of the matrixes
m = zeros(degree+1);
l = zeros(degree+1);
M = zeros((degree+1)+n*degree-1);
L = zeros((degree+1)+n*degree-1);
D = zeros((degree+1)+n*degree-1);
a = 0;
b = degree*h;
for j = 1:degree+1
    for i = 1:degree+1
        m(i,j) = diff(polyval(polyint(conv(phi(j,:),phi(i,:))),[a,b]));
        l(i,j) = diff(polyval(polyint(conv(phiPrim(j,:),phi(i,:))),[a,b]));
        d(i,j) = diff(polyval(polyint(conv(phiPrim(j,:),phiPrim(i,:))),[a,b]));
    end
end

for i = 1:degree:length(M)-degree % Assemble FEM matrixes
    M(i:degree+i,i:degree+i) = M(i:degree+i,i:degree+i) + m;
    L(i:degree+i,i:degree+i) = L(i:degree+i,i:degree+i) + l;
    D(i:degree+i,i:degree+i) = D(i:degree+i,i:degree+i) + d;
end


% to fix last addition of local 
M(1,1) = M(1,1) + m(end,end);
L(1,1) = L(1,1) + l(end,end);
D(1,1) = D(1,1) + d(end,end);

M(end-degree+1:end,end-degree+1:end) = M(end-degree+1:end,end-degree+1:end) + m(1:degree,1:degree);
L(end-degree+1:end,end-degree+1:end) = L(end-degree+1:end,end-degree+1:end) + l(1:degree,1:degree);
D(end-degree+1:end,end-degree+1:end) = d(end-degree+1:end,end-degree+1:end) + d(1:degree,1:degree);

M(1,end-degree+1:end) = M(1,end-degree+1:end) + m(end,1:degree);
L(1,end-degree+1:end) = L(1,end-degree+1:end) + l(end,1:degree);
D(1,end-degree+1:end) = D(1,end-degree+1:end) + d(end,1:degree);

M(end-degree+1:end,1) = M(end-degree+1:end,1) + m(1:degree,end);
L(end-degree+1:end,1) = L(end-degree+1:end,1) + l(1:degree,end);
D(end-degree+1:end,1) = D(end-degree+1:end,1) + d(1:degree,end);



% print matrixes
M % mass matrix
L % phiPrim * phi matrix
D % phiPRim * phiPrim matrix (stabilization)
