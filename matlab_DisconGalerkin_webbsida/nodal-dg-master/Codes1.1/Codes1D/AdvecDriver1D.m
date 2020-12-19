function AdvecDriver1D(degrees,bounds,intervals, alpha)
% degrees : the Order of polymomials used for approximation
% bounds : the left and the right boundaries t ex [0 1]
% intervals: the number of elements
% alpha:  numerical flux (stable: 0<=alpha<=1)
% alpha = 1; central flux alpha =0; % upwind flux


% Driver script for solving the 1D advection equations
Globals1D;

% Order of polymomials used for approximation
%N = 8; %%%%%%%%%%%%%%%%changed 
N=degrees;

% Generate simple mesh
% [Nv, VX, K, EToV] = MeshGen1D(0.0,2.0,10); %%%%%%%%%%%%%%%%changed
% [Nv, VX, K, EToV] = MeshGen1D(0.0,1.0,10);
[Nv, VX, K, EToV] = MeshGen1D(bounds(1),bounds(2),intervals);
% Nv % nr of points
% VX node coordinates
% K nr of elements
% EToV a matrix containg the nr of the nodes that a specific element has ie
% x1 to x2 is the first rad (the first element)

% Initialize solver and construct grid and metric
StartUp1D;
% Set initial conditions
% u = sin(x);  %%%%%%%%%%%%%%%%changed
%% analytic solution
u_0 = 1; % amplitude
k = 2*pi; % wave frequency
u = real(u_0*exp(1i*k*(x)));
%%
% advection speed  %%%%%%%%%%%%%%%%changed
%a = 2*pi;
a=1;
% numerical flux (stable: 0<=alpha<=1)
%alpha = 1; % central flux
% alpha = 0; % upwind flux
% Solve Problem
FinalTime = 10;
[u] = Advec1D(u,FinalTime,a,alpha);
end