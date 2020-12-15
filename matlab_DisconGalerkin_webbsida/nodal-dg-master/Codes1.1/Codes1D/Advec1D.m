function [u] = Advec1D(u, FinalTime, a, alpha)
% Purpose : Integrate 1D advection until FinalTime starting with
% initial the condition, u
Globals1D;

% CFL constant
CFL=0.75;

time = 0;
% Runge-Kutta residual storage
resu = zeros(Np,K);
% compute time step size
dxmin = min(abs(x(1,:)-x(2,:)));
dt = CFL*dxmin/a;
Nsteps = ceil(FinalTime/dt); dt = FinalTime/Nsteps;
% outer time step loop
for tstep=1:Nsteps
    for INTRK = 1:5
        timelocal = time + rk4c(INTRK)*dt;
        [rhsu] = AdvecRHS1D(u, timelocal, a, alpha);
        resu = rk4a(INTRK)*resu + dt*rhsu;
        u = u+rk4b(INTRK)*resu;
    end;
% Increment time
time = time+dt;
% Plot solution
plot(x,u); drawnow;
end;
return