% Driver script for solving the 1D advection equations
Globals1D;
% Order of polymomials used for approximation
N = 8;
% Generate simple mesh
[Nv, VX, K, EToV] = MeshGen1D(0.0,2.0,10);
% Nv % nr of points
% VX node coordinates
% K nr of elements
% EToV a matrix containg the nr of the nodes that a specific element has ie
% x1 to x2 is the first rad (the first element)

% Initialize solver and construct grid and metric
StartUp1D;
% Set initial conditions
u = sin(x);
% advection speed
a = 2*pi;
% numerical flux (stable: 0<=alpha<=1)
%alpha = 1; % central flux
alpha = 0; % upwind flux
% Solve Problem
FinalTime = 10;
[u] = Advec1D(u,FinalTime,a,alpha);