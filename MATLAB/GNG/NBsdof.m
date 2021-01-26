function [ x, xd, xdd ] = NBsdof( f, m, c, k, dt, x0, xd0, xdd0 )
%NBsdof Constant average acceleration Newmark-Beta numerical integration
%algorithm for single degree of freedom system

x = (f + m*(4*x0/dt^2 + 4*xd0/dt + xdd0) + c*(2*x0/dt + xd0)) / (4*m/dt^2 + 2*c/dt + k);

xd = -xd0 + 2*(x - x0)/dt;

xdd = 4*(x - x0 - xd0*dt)/dt^2 - xdd0;

end

