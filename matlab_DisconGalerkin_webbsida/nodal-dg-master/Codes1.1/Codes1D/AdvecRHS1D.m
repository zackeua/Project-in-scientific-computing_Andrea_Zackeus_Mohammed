function [rhsu] = AdvecRHS1D(u,time, a,alpha)
% function [rhsu] = AdvecRHS1D(u,time)
% Purpose : Evaluate RHS flux in 1D advection
Globals1D;
% form flux differences at faces
df = zeros(Nfp*Nfaces,K);
df(:) = 0.5*a*(u(vmapM)-u(vmapP)).*(nx(:)-(1-alpha));
% impose boundary condition at x=0
% uin = -sin(a*time);  %%%%%%%%%%%%%%%%changed
% i let BC =exact sol
%% analytic solution
u_0 = 1; % amplitude
k = 2*pi; % wave frequency
uin = real(u_0*exp(1i*k*(-time)));

%
df(mapI) = 0.5*a*(u(vmapI)-uin).*(nx(mapI)-(1-alpha));
df(mapO) = 0;
% compute right hand sides of the semi-discrete PDE
rhsu = -a*rx.*(Dr*u) + LIFT*(Fscale.*(df));
return