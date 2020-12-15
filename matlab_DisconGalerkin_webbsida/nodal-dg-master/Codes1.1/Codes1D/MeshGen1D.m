function [Nv, VX, K, EToV] = MeshGen1D(xmin,xmax,K)
%fattar bra
% copied from ServiceRoutines by me
% function [Nv, VX, K, EToV] = MeshGen1D(xmin,xmax,K)
% Purpose  : Generate simple equidistant grid with K elements

Nv = K+1; % nr of points

% Generate node coordinates
VX = (1:Nv);
for i = 1:Nv
  VX(i) = (xmax-xmin)*(i-1)/(Nv-1) + xmin; % gives 0, 0.2, 0.4...
end

% read element to node connectivity
EToV = zeros(K, 2);
for k = 1:K
  EToV(k,1) = k; EToV(k,2) = k+1;
  % EToV gives EToV =

%      1     2
%      2     3
%      3     4
%      4     5
%      5     6
%      6     7
%      7     8
%      8     9
%      9    10
%     10    11 ie x1 to x2 the first rad: makes the first element and so on

end
return
