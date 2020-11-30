function [M,L,K,X] = MatrixAssembler(degree,n,mode)
% degree is degree of polynomials,
% n is the number of intervals
% mode is the type of matrix assembly
% 1 = exact integration, 2=gauss lobatto integration

local_M = zeros(degree+1);
local_L = zeros(degree+1);
local_K = zeros(degree+1);

M = zeros(degree*(n+1));
L = zeros(degree*(n+1));
K = zeros(degree*(n+1));

phiPrim = zeros(degree+1,degree);

if mode==1
    X = 0:1/n:1;
    phi = coeff2(degree,X(1):X(2)/degree:X(2));
end

if mode==2
    [X,W]= legendre_gauss_lobatto(degree+1);
    right = 1/(n+1);
    left = 0;
    X = (right-left)/2 * X +(right + left)/2;
    X = flip(X);
    W  = W*(right-left)/2;
    phi = coeff2(degree,X);
end

if mode==3
    [X,W]= legendre_gauss_lobatto(degree+1);
    right = 1/(n+1);
    left = 0;
    X = (right-left)/2 * X +(right + left)/2;
    X = flip(X);
    W  = W*(right-left)/2;
    phi = coeff2(degree,X);
end



for i = 1:degree+1
    phiPrim(i,:) = polyder(phi(i,:));
end


if mode==1
    a = X(1);
    b = X(2);

    h = 1/degree/(n+1);
    X = zeros(degree*(n+1),1);
    for i = 1:degree*(n+1)
        X(i) = h*(i-1);
    end
    X = linspace(0,1,degree*(n+1)+1);
    X = X(1:end-1)';

    for j = 1:degree+1
        for i = 1:degree+1
            local_M(i,j) = diff(polyval(polyint(conv(phi(j,:),phi(i,:))),[a,b]));
            local_L(i,j) = diff(polyval(polyint(conv(phiPrim(j,:),phi(i,:))),[a,b]));
            local_K(i,j) = diff(polyval(polyint(conv(phiPrim(j,:),phiPrim(i,:))),[a,b]));
        end
    end
end

if mode==2
    for j = 1:degree+1
        for i = 1:degree+1
            f1 = conv(phi(j,:),phi(i,:));
            f1x = polyval(f1,X);
            local_M(i,j) = W'*f1x;
            
            f2 = conv(phiPrim(j,:),phi(i,:));
            f2x = polyval(f2,X);
            local_L(i,j) = W'*f2x;
            
            f3 = conv(phiPrim(j,:),phiPrim(i,:));
            f3x = polyval(f3,X);
            local_K(i,j) = W'*f3x;
            
        end
    end
    Xout = X;
    for i = 1:n
        Xout = [Xout; Xout(end) + X(2:end)];
    end
    X = Xout(1:end-1);
end


if mode==3
    a = X(1);
    b = X(end);
    
    for j = 1:degree+1
        for i = 1:degree+1
            local_M(i,j) = diff(polyval(polyint(conv(phi(j,:),phi(i,:))),[a,b]));
            local_L(i,j) = diff(polyval(polyint(conv(phiPrim(j,:),phi(i,:))),[a,b]));
            local_K(i,j) = diff(polyval(polyint(conv(phiPrim(j,:),phiPrim(i,:))),[a,b]));
        end
    end
    Xout = X;
    for i = 1:n
        Xout = [Xout; Xout(end) + X(2:end)];
    end
    X = Xout(1:end-1);
end




for i = 1:degree:degree*n
    M(i:i+degree,i:i+degree) = M(i:i+degree,i:i+degree) + local_M;
    L(i:i+degree,i:i+degree) = L(i:i+degree,i:i+degree) + local_L;
    K(i:i+degree,i:i+degree) = K(i:i+degree,i:i+degree) + local_K;
end

% to fix last addition of local 
M(1,1) = M(1,1) + local_M(end,end);
L(1,1) = L(1,1) + local_L(end,end);
K(1,1) = K(1,1) + local_K(end,end);

M(end-degree+1:end,end-degree+1:end) = M(end-degree+1:end,end-degree+1:end) + local_M(1:degree,1:degree);
L(end-degree+1:end,end-degree+1:end) = L(end-degree+1:end,end-degree+1:end) + local_L(1:degree,1:degree);
K(end-degree+1:end,end-degree+1:end) = K(end-degree+1:end,end-degree+1:end) + local_K(1:degree,1:degree);

M(1,end-degree+1:end) = M(1,end-degree+1:end) + local_M(end,1:degree);
L(1,end-degree+1:end) = L(1,end-degree+1:end) + local_L(end,1:degree);
K(1,end-degree+1:end) = K(1,end-degree+1:end) + local_K(end,1:degree);

M(end-degree+1:end,1) = M(end-degree+1:end,1) + local_M(1:degree,end);
L(end-degree+1:end,1) = L(end-degree+1:end,1) + local_L(1:degree,end);
K(end-degree+1:end,1) = K(end-degree+1:end,1) + local_K(1:degree,end);





end