function[x1, xd1, xdd0] = dynamic1(p0, dp, c0, k0, x0, xd0, M, dt)
%NUMERICAL INTEGRATION
%   Newmark-beta incremental integration algorithm

%current damping and spring forces
fd0 = c0*xd0;
fs0 = k0*x0;

%initial acceleration
xdd0 = ((p0 - fd0 - fs0)' / M)';

%effective stiffness
kc = k0 + 2*c0/dt + 4*M/(dt^2);

%effective load increment
dpc = dp + M*(4*xd0/dt + 2*xdd0) + 2*c0*xd0;

%displacement and velocity increments
dx = kc \ dpc;
dxd = 2*dx/dt - 2*xd0;

%final velocity and displacement
xd1 = xd0 + dxd;
x1 = x0 + dx;

end

