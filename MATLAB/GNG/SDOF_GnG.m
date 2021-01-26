% SDOF GnG Hysteresis Model
% jco118

close all; clear; clc; %reset
set(0,'units','centimeters','defaultfigureposition',[60 60 840 660]); %new default figure

m = 1e7;
c = 2e8;
k = 4e6;%4e6;%3e7;%400;

load('earthquake signals\Standard Set\t_0_001s\EQ_records_Low.mat');
endgame = 10000;
equake = Accel_matrix(1:endgame-1,2); % Coyote Lake, 1979, Gilroy Array 2 (la42)
acc = [0; equake/1];
t = Time_vector(1:endgame);
dt = t(1);

% dt = 0.01;
% t = (0:dt:250)';
% acc = -0.1*sin(0.7*t);
f = -m*acc;
N = length(t);
fGnG = zeros(N,1);

% Response Preallocation
% Nominal structure
x = zeros(N,1); % displacement
xd = zeros(N,1); % velocity
xdd = zeros(N,1); % acceleration
% Nominal structure with devices added
xGnG = zeros(N,1); % displacement
xdGnG = zeros(N,1); % velocity
xddGnG = zeros(N,1); % acceleration

% INITIALISE GNG PROPERTIES
GnG.E = 200e9; % elastic modulus
GnG.d = 0.20; % fuse diameter (0.0116)
GnG.l = 0.7; % system length
GnG.fl = 0.3; % fuse length
GnG.sigmay = 393e6; % yield stress (fuse)
GnG.K2 = 6e4; % post-yield stiffness (system)
GnG.Beta = 8; % hysteresis divergence coefficient (working range of up to 140)
GnG.pitch = 20e-3; % rack pitch
GnG.eab = 0.28; % fuse elongation at break
GnG.teeth = 300; % available teeth on rack

GnG.fusebreak = GnG.fl*(1+GnG.eab); % fuse length at fracture
GnG.A = (pi/4)*GnG.d^2; % cross-section area (fuse)
GnG.K1 = GnG.E*GnG.A/GnG.l; % initial stiffness (system)
GnG.Fy = GnG.sigmay*GnG.A; % yield force (fuse)
GnG.xy = GnG.Fy/GnG.K1; % yield displacement
GnG.Fr = 0; % reset force (hysteresis)
GnG.xr = 0; % reset displacement (hysteresis)
GnG.Fp = 0; % previous force (hysteresis)
GnG.xp = 0; % previous displacement (hysteresis)
GnG.x0 = 0; % engagement displacement (hysteresis)
GnG.ratchet = []; % array of ratcheting time steps
GnG.sx_p = 1; % sign of previous displacement step (defaulting to positive)

GnG.F = zeros(N,1); % force
GnG.fuse = ones(N,1)*GnG.fl; % fuse length
GnG.dx = zeros(N,1); % displacement step

device1 = GnG;
device2 = GnG;

tic;

for ii = 2:N
    
    % Response of nominal structure
    
    [ x(ii), xd(ii), xdd(ii) ] = ...
        NBsdof( f(ii), m, c, k, dt, x(ii-1), xd(ii-1), xdd(ii-1) );
    
    % Response of nominal structure with devices added
    
    sxd = sign(xdGnG(ii-1));
    fGnG(ii) = f(ii) - sxd*device1.F(ii-1) + sxd*device2.F(ii-1);
    
    [ xGnG(ii), xdGnG(ii), xddGnG(ii) ] = ...
        NBsdof( fGnG(ii), m, c, k, dt, xGnG(ii-1), xdGnG(ii-1), xddGnG(ii-1) );
    
    % Response of GnG devices
    
    [ device1 ] = GnGfun( xGnG(ii), device1, ii );
    [ device2 ] = GnGfun( -1*xGnG(ii), device2, ii );
    
end

toc;

figure, plot(t,[x xd xdd])
legend('Displacement','Velocity','Acceleration')
title('Nominal Structure Response'), xlabel('Time (s)')
ylabel('Displacement (m), Velocity (m/s), Acceleration (m/s^2)')

figure, plot(t,[xGnG xdGnG xddGnG])
legend('Displacement','Velocity','Acceleration')
title('Instrumented Structure Response'), xlabel('Time (s)')
ylabel('Displacement (m), Velocity (m/s), Acceleration (m/s^2)')

figure, plot(xGnG,[device1.F device2.F])
legend('Device 1','Device 2')
title('GnG Device Hysteresis'), xlabel('Displacement'), ylabel('Force')

% for jj = 1:N
%     
%     h = figure(1);
%     plot(x(1:jj),device2.F(1:jj),'g-',x(jj),device2.F(jj),'r+')
%     
% end

% figure, plot(device1.fuse-320e-3,device1.F)
% figure, plot(device2.fuse-320e-3,device2.F)
% 
% figure, plot(x,device1.F,device1.fuse-320e-3,device1.F)
% figure, plot(x,device2.F,device2.fuse-320e-3,device2.F)
% figure, plot(device1.fuse-320e-3,device1.F,device2.fuse-320e-3,device2.F)

figure, plot(t,[x xd xdd],t,[xGnG xdGnG xddGnG])
legend('Displacement (nominal)','Velocity (nominal)','Acceleration (nominal)',...
    'Displacement (instrumented)','Velocity (instrumented)','Acceleration (instrumented)')
title('Nominal and Instrumented Structure Response Comparison'), xlabel('Time (s)')
ylabel('Displacement (m), Velocity (m/s), Acceleration (m/s^2)')

figure, plot(t,[device1.F device2.F f fGnG])
title('Force Comparison'), xlabel('Time (s)'), ylabel('Force (N)')
legend('GnG 1','GnG 2','EQ Load','Structure')

max_values_nominal = max([abs(x) abs(xd) abs(xdd)]);
display(max_values_nominal);
max_values_instrumented = max([abs(xGnG) abs(xdGnG) abs(xddGnG)]);
display(max_values_instrumented);
response_percent_w_GnG = max_values_instrumented./max_values_nominal*100;
display(response_percent_w_GnG);
