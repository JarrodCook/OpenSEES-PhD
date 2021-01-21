%% Load OpenSees data from leaning column
% clear previous work
close all; clear; clc;

%% Schematic

%                          N2                           ^
%                          |                            |
%                          |                            |
%                          |                            |
%                          |                            |
%                          |                            |
%                          |                            |
%                          E1                           |
%                          |                            |
%                          |                            | H
%                          |                            |
%                          |                            |
%                          |                            |
%                          |                            |
%      ---E6---N3----E2----N1----E3----N4---E7---       |
%      |       |           |           |        |       |
%      |       E4          E8          E5       |       |
%      |       |           |           |        |       |
%      --------N5          N7          N6--------       v
%
%              <-----------><----------->
%                    b            b
% Nodes
% N1 = frame base node
% N2 = frame roof node
% N3, N4 = rocking edge nodes
% N5, N6, N7 = fixed ground nodes
% Elements
% E1 = frame column element
% E2, E3 = frame beam elements
% E4, E5 = zero length rocking and device elements
% E6, E7 = zero length horizontal support elements
% E8 = zero length PT element

%% run .tcl file in OpenSees

% AccelRamp = linspace(0,20,20e3);
% save AccelRamp.txt AccelRamp -ascii

!OpenSeesGNG.exe "run_thisIsIt".tcl
% % % !OpenSeesGNG.exe "run_frameOnly".tcl

% return

%% load in the data

filedir = 'Analysis_Results/thisIsIt/';

% test details
testDetails     = load(strcat(filedir,'testDetails.txt'));

EQnum = testDetails(1);
B = 2.5;
AR = testDetails(5);
H = 2*B*AR;
HM = 2*H/3;
W = 1e5*9.81;
MAll = 21e5;

% node recorders
DispData        = load(strcat(filedir,'Disps.txt'));
AccelData       = load(strcat(filedir,'Accs.txt'));
ReactionData    = load(strcat(filedir,'Reactions.txt'));

% element recorders
ForceData       = load(strcat(filedir,'Forces.txt'));
DefData         = load(strcat(filedir,'Defs.txt'));
DemandData      = load(strcat(filedir,'Demand.txt'));
RatchetData     = load(strcat(filedir,'Ratchet.txt'));

%% extract data to vectors

% time and displacements
time    = DispData(:,1);
dispx   = DispData(:,2:3:end);
dispy   = DispData(:,3:3:end);
dispz   = DispData(:,4:3:end);

% accelerations
accelx  = AccelData(:,1:3:end);
accely  = AccelData(:,2:3:end);
accelz  = AccelData(:,3:3:end);

% reactions
reactx  = ReactionData(:,1:3:end);
reacty  = ReactionData(:,2:3:end);
reactz  = ReactionData(:,3:3:end);

% forces (Careful - not applicable to all elements!!)
forcex  = ForceData(:,1:6:end);
forcey  = ForceData(:,2:6:end);
forcez  = ForceData(:,3:6:end);

% deformations (Careful - not applicable to all elements!!)
defx    = DefData(:,1:6:end);
defy    = DefData(:,2:6:end);

% demand
demand = DemandData; % LHS = index 1, RHS = index 2
ratchet = RatchetData; % LHS = index 1, RHS = index 2

%% figures
% [xstart ystart xdim ydim] - (reg. size ~ 560x420)
% Uplift at rocking edges
% figure('Position',[50 550 560 420])
figure('Position',[50 200 1800 800])
subplot(2,5,3)
plot(time,[dispy(:,3),dispy(:,4)])
grid on
xlabel('Time')
ylabel('Displacement')
title('RE vertical disp.')
% set(gca,'fontsize',16)
legend('Left foot','Right foot')
% Force-displacement (roof force)
% figure('Position',[650 550 560 420])
subplot(2,5,5)
plot(dispx(:,2)-dispx(:,1),reactx(:,9)+reactx(:,2))
grid on
xlabel('Lateral deflection')
ylabel('Roof force')
title('At mass height')
% set(gca,'fontsize',16)
% Force-displacement (base force)
% figure('Position',[1250 450 560 420])
subplot(2,5,4)
plot(dispx(:,2)-dispx(:,1),[-reactx(:,5),-reactx(:,6)])
grid on
xlabel('Lateral deflection')
ylabel('Base reaction force')
title('Force-displacement at mass')
% set(gca,'fontsize',16)
legend('Left foot','Right foot','Location','Best')

% return

% Rocking edge force-displacement plots (vertical and horizontal)
% figure('Position',[50 50 700 500])
% suptitle('Rocking edge force-displacements results')
% Vertical - left
ax1 = subplot(2,5,1); % top left
plot(dispy(:,3),-reacty(:,5),'bo-','MarkerSize',2)
grid on
title('LHS GNG - vertical')
xlabel('Displacement')
ylabel('Reaction force')
% axis([min(min(dispy(:,3:4))) max(max(dispy(:,3:4))) min(min(-reacty(:,5:6))) max(max(-reacty(:,5:6)))])
% Vertical - right
ax2 = subplot(2,5,2); % top right
plot(dispy(:,4),-reacty(:,6),'ro-','MarkerSize',2)
grid on
title('RHS GNG - vertical')
xlabel('Displacement')
ylabel('Reaction force')
% linkaxes([ax1,ax2],'xy')
% Horizontal - left
ax3 = subplot(2,5,6); % bottom left
plot(dispx(:,3),-reactx(:,5),'bo-','MarkerSize',2)
grid on
title('LHS GNG - horizontal')
xlabel('Displacement')
ylabel('Reaction force')
% Horizontal - right
ax4 = subplot(2,5,7); % bottom right
plot(dispx(:,4),-reactx(:,6),'ro-','MarkerSize',2)
grid on
title('RHS GNG - horizontal')
xlabel('Displacement')
ylabel('Reaction force')

% figure('Position',[1250 50 560 420])
% plot(time,DamageData)
% grid on
% xlabel('time')
% ylabel('damage')
% set(gca,'fontsize',16)
% legend('Left device','Right device','Location','NorthWest')

% % dissipater elongation calcs! ***(assumes no load reversal mid rock?)***
% dyield = 6e-5; % yield displacement
% pksleft = findpeaks(dispy(:,3)); % peak displacements of LHS device
% pksright = findpeaks(dispy(:,4)); % peak displacements of RHS device
% bigpkslocleft = pksleft > dyield; % indices of peaks beyond yield - LHS
% bigpkslocright = pksright > dyield; % indices of peaks beyond yield - RHS
% bigpksleft = pksleft(bigpkslocleft); % peak values beyond yield - LHS
% bigpksright = pksright(bigpkslocright); % peak values beyond yield - RHS
% elongleft = bigpksleft - dyield; % steps of elongation - LHS
% elongright = bigpksright - dyield; % steps of elongation - RHS
% totalelongleft = sum(elongleft); % total elongation - LHS
% totalelongright = sum(elongright); % total elongation - RHS
% fprintf('Elongation of LHS dissipater = %5.2f mm\n',totalelongleft*1000)
% fprintf('Elongation of RHS dissipater = %5.2f mm\n',totalelongright*1000)

% figure
subplot(2,5,8)
plot(dispy(:,1),-reacty(:,7),'bo-','MarkerSize',2)
grid on
xlabel('Displacement')
ylabel('Force')
title('Post-tensioning')

subplot(2,5,9)
plot(time,dispx(:,10))
grid on
xlabel('Time')
ylabel('Displacement')
title('Roof displacement')

subplot(2,5,10)
plot(time,demand)
grid on
xlabel('Time')
ylabel('Demand (m)')
title('Inelastic dissipater demand')
legend('GNG LHS','GNG RHS','Location','Best')

% figure
% plot(dispx(:,10),[-reactx(:,5),-reactx(:,6)])
% grid on
% xlabel('Lateral deflection')
% ylabel('Base reaction force')
% title('System force-displacement')
% % set(gca,'fontsize',16)
% legend('Left foot','Right foot','Location','Best')

% Print analysis details
fprintf('Test details:\n')
fprintf('EQ: %i, Period: %.2f s, Scale factor: %.4f, Gup: %.4f g\n',testDetails(1:4))
fprintf('Aspect ratio: %i, Pitch: %.3f m, PT: %i %%, R: %i\n',testDetails(5:8))
fprintf('Yield force: %.4f N, OK: %i, Eigenvalue: %.4f\n',testDetails(9:11))

% Energy dissipation
yieldForce = testDetails(9);
energyDissLHS = yieldForce*demand(end,1);
energyDissRHS = yieldForce*demand(end,2);

fprintf('\nEnergy dissipation (LHS) = %.2f Nm\n',energyDissLHS)
fprintf('Energy dissipation (RHS) = %.2f Nm\n',energyDissRHS)

% Maxmium deflections at base
LHSdefy = min(dispy(:,3))*1000;
RHSdefy = min(dispy(:,4))*1000;
Cdefy = min(dispy(:,1))*1000;
LHSdefx = max(abs(dispx(:,3)))*1000;
RHSdefx = max(abs(dispx(:,4)))*1000;

fprintf('\nMaximum downward deflection at base:\n')
fprintf('LHS: %.3f mm, Centre: %.3f mm, RHS: %.3f mm\n',LHSdefy,Cdefy,RHSdefy)
fprintf('\nMaximum horizontal deflection at base:\n')
fprintf('LHS: %.3f mm, RHS: %.3f mm\n',LHSdefx,RHSdefx)

% figure
% plot(dispx(:,2),reactx(:,2))
% grid on

% figure
% plot(time,dispx(:,2))
% grid on
% ylabel('Disp. at mass')
% figure
% plot(time,reactx(:,2))
% grid on
% ylabel('Force at mass')
% % figure
% % plot(dispx(:,2),reactx(:,2))
% % figure
% % plot(time,dispx(:,9))
% % grid on
% % ylabel('Disp. at trib. mass')
% figure
% plot(time,reactx(:,9))
% grid on
% ylabel('Force at trib. mass')
% % figure
% % plot(dispx(:,9),reactx(:,9))
% % grid on

% Tn = 0.5;
% MAll = 21e5;
% Kinitial = (2*pi/Tn)^2*MAll;
% disp('Kinitial:')
% disp(Kinitial)
% 
% feauxK = max(-reactx(:,6))/max(dispx(:,10));
% disp('Feaux K: (N/A after uplift)')
% disp(feauxK)
% 
% KatLC = max(reactx(:,9))/max(dispx(:,9));
% disp('K at leaning column: (unreliable)')
% disp(KatLC)

% % Hysteresis comparison
% figure
% plot(dispx(:,2),reactx(:,9)+reactx(:,2),'g-')
% hold on
% plot(dispx(:,2),[-reactx(:,5) -reactx(:,6)])
% grid on
% xlabel('Displacement (m)')
% ylabel('Force (N)')
% legend('At Mass','Left Foot','Right Foot','Location','Best')
% title('Hysteresis comparison')

% figure
% plot(time,[dispy(:,3) dispy(:,4) dispy(:,1)])
% grid on
% legend('Left','Right','Centre')
% xlabel('Time (s)')
% ylabel('Vertical displacement (m)')

% Convert base rotation to deflection
% B = 2.5;
% AR = 1;
% H = 2*B*AR;
% HM = H*2/3;

defFromRot = dispz(:,1)*HM;

% figure
% plot(time,[dispx(:,2) -defFromRot])
% grid on


% figure
% plot(dispx(:,2)-dispx(:,1)+defFromRot,[-reactx(:,5) -reactx(:,6)])
% grid on
% xlabel('Displacement minus that from base rotation (m)')
% ylabel('Force at rocking edge (N)')
% title('Column stiffness')

Ktrue = max(abs([-reactx(:,5);-reactx(:,6)])) / max(abs(dispx(:,2)+defFromRot-dispx(:,1)));
fprintf('\nKtrue: %.4E\n',Ktrue)

maxDrift = max(abs(dispx(:,10)))/H;
fprintf('Maximum drift (old) = %.4f %%\n',maxDrift*100)
maxDriftTrue = max(abs(dispx(:,10)-dispx(:,1)))/H;
fprintf('Maximum drift (true) = %.4f %%\n',maxDriftTrue*100)

% close all

% figure
% plot(dispx(:,2)-dispx(:,1)+defFromRot,-reactx(:,5))
% grid on

% return

% Hysteresis comparison
% figure
% plot(dispx(:,2)-dispx(:,1),reactx(:,9)+reactx(:,2),'g-')
% grid on
% hold on
% plot(dispx(:,2)-dispx(:,1),[-reactx(:,5) -reactx(:,6)])
% legend('Force at mass','Left edge','Right edge','Location','Best')
% title('Finally, the real hysteresis plot! (at mass height)')
% xlabel('Displacement (m)')
% ylabel('Force at rocking edge (N)')

KsLHS = ( reactx(2:end,5) - reactx(1:end-1,5) ) ./ ...
    ( (dispx(2:end,2)-dispx(2:end,1)) - (dispx(1:end-1,2)-dispx(1:end-1,1)) );
KsRHS = ( reactx(2:end,6) - reactx(1:end-1,6) ) ./ ...
    ( (dispx(2:end,2)-dispx(2:end,1)) - (dispx(1:end-1,2)-dispx(1:end-1,1)) );

% Deflection at top of frame
figure
plot(dispx(:,10)-dispx(:,1),[-reactx(:,5) -reactx(:,6)])
grid on
xlabel('Displacement at Node 10 (m)')
ylabel('Rocking base node horizontal force (MN)')
% title('Deflection at top of frame')

figure
plot(dispy(:,3),-reacty(:,5),'bo-','MarkerSize',2)
grid on
title('LHS GNG - vertical')
xlabel('Displacement')
ylabel('Reaction force')

figure
plot(dispy(:,4),-reacty(:,6),'ro-','MarkerSize',2)
grid on
title('RHS GNG - vertical')
xlabel('Displacement')
ylabel('Reaction force')

% Hysteresis
figure
plot((dispx(:,10)-dispx(:,1))/H*100,[-reactx(:,5) -reactx(:,6)]/1e6)
grid on
axis([-1.5 1.5 -8 8])
xlabel('Drift at Node 10 (%)')
ylabel('Rocking base node horizontal force (MN)')
legend('LHS rocking edge','RHS rocking edge','Location','Best')
xt = get(gca, 'XTick');
set(gca, 'FontSize', 18)

% Deflection at top of frame
figure
plot(dispx(:,10)-dispx(:,1),[-reactx(:,5) -reactx(:,6)]/1e6)
grid on
axis([-0.25 0.25 -8 8])
xlabel('Displacement at Node 10 (m)')
ylabel('Rocking base node horizontal force (MN)')
legend('LHS rocking edge','RHS rocking edge','Location','Best')
% title('Deflection at top of frame')
xt = get(gca, 'XTick');
set(gca, 'FontSize', 18)

figure
plot(time,[dispy(:,3),dispy(:,4)])
grid on
xlabel('Time')
ylabel('Displacement')
title('RE vertical disp.')
% set(gca,'fontsize',16)
legend('Left foot','Right foot')

return

%% Figures towards comment S-15!!

% RESPONSE

% Drift at top of frame
figure
plot(time,(dispx(:,10)-dispx(:,1))/H*100)
grid on
xlabel('Time (s)')
ylabel('Drift at node 10 (%)')
xt = get(gca, 'XTick');
set(gca, 'FontSize', 22)
fileLoc = 'FiguresS15/Drift';
saveas(gcf,fileLoc,'meta')
saveas(gcf,fileLoc,'fig')

% Uplift at rocking edges
figure
plot(time,[dispy(:,3),dispy(:,4)]*1000)
grid on
ylim([-1 80])
xlabel('Time (s)')
ylabel('Uplift (mm)')
legend('LHS rocking edge','RHS rocking edge')
xt = get(gca, 'XTick');
set(gca, 'FontSize', 22)
fileLoc = 'FiguresS15/Uplift';
saveas(gcf,fileLoc,'meta')
saveas(gcf,fileLoc,'fig')

% Hysteresis
figure
plot((dispx(:,10)-dispx(:,1))/H*100,[-reactx(:,5) -reactx(:,6)]/1e6)
grid on
axis([-2 2 -8 8])
xlabel('Drift at Node 10 (%)')
ylabel('Rocking base node horizontal force (MN)')
legend('LHS rocking edge','RHS rocking edge','Location','Best')
xt = get(gca, 'XTick');
set(gca, 'FontSize', 22)
fileLoc = 'FiguresS15/FD';
saveas(gcf,fileLoc,'meta')
saveas(gcf,fileLoc,'fig')

% GROUND MOTION

scaleFactors = load('scaleFactors.txt');
AllEQs = load('earthquake signals/Standard Set/t_0_001s/EQ_records_All.mat');
timeGM = [0 AllEQs.Time_vector];
Accels = AllEQs.Accel_matrix;
NGM = size(Accels,2);
Accels = [zeros(1,NGM); Accels];
Ts = (0.05:0.05:5);
NT = length(Ts);
GMno = 1;
thePeriod = testDetails(2);
TnRef = find(Ts == thePeriod);
TnRef = 14;

% FANCY COMMON TIME AXIS PLOT

figure('Position',[50 200 1200 800])
ax1 = subplot(3,1,1);
plot(timeGM,scaleFactors(TnRef,GMno)*Accels(:,GMno))
grid on
ylim([-4 4])
xlim([0 50])
ylabel({'Acceleration','(ms^-^2)'})
xt = get(gca, 'XTick');
set(gca, 'FontSize', 22)
ax2 = subplot(3,1,2);
plot(time,(dispx(:,10)-dispx(:,1))/H*100)
grid on
ylim([-2 2])
xlim([0 50])
ylabel({'Drift at node 10','(%)'})
xt = get(gca, 'XTick');
set(gca, 'FontSize', 22)
ax3 = subplot(3,1,3);
plot(time,[dispy(:,3),dispy(:,4)]*1000)
grid on
ylim([-inf 80])
xlim([0 50])
ylabel({'Uplift','(mm)'})
xlabel('Time (s)')
legend('LHS rocking edge','RHS rocking edge')
xt = get(gca, 'XTick');
set(gca, 'FontSize', 22)
linkaxes([ax1,ax2,ax3],'x')
fileLoc = 'FiguresS15/GM01DriftUplift';
saveas(gcf,fileLoc,'meta')
saveas(gcf,fileLoc,'fig')

% legend('LHS rocking edge (No GNG)','RHS rocking edge (No GNG)',...'LHS rocking edge (GNG)','LHS rocking edge (GNG)','Location','EastOutside')

%% Reconstructing the hysteresis flag (towards amendment S-20!)

g = 9.81;
R = 4;
PTper = 100;
Tn = 0.7;
wn = 2*pi/Tn;
Kinitial = wn*wn*MAll;
WAll = MAll*g;
myAlpha = 0.8;
myGamma = 0.1;

maxD = max(abs(dispx(:,10)));
maxF = max(abs([reactx(:,5);reactx(:,6)]));

AA = load('Gups.txt');
Gup = AA(14);

Mover = Gup*g*MAll*HM/R;
MW = W*B;
MED = Mover*(100-PTper)/100;
MPT = Mover - MW - MED;

y1 = (MW + MPT)/HM;
x1 = y1/Kinitial;

y2 = Mover/HM;
x2 = x1 + (MED/HM)/(Kinitial*myAlpha);

y3 = maxF;
x3 = maxD;

y4 = y3 - (y2-y1);
x4 = x3 - (x2-x1);

flagx = [x1 x2 x3 x4 x1];
flagy = [y1 y2 y3 y4 y1];

% The plot

figure
plot(dispx(:,10)-dispx(:,1),-reactx(:,5)-reactx(:,6))
grid on
% axis([-1.5 1.5 -8 8])
xlabel('Displacement at Node 10 (m)')
ylabel('Rocking base node horizontal force (N)')
legend('LHS rocking edge','RHS rocking edge','Location','Best')
xt = get(gca, 'XTick');
set(gca, 'FontSize', 18)

hold on
plot(flagx,flagy,'m--')

linearD = maxF*maxD/2;

widthD = sqrt((x3-x2)^2+(y3-y2)^2);
heightD = sqrt((x2-x1)^2+(y2-y1)^2);
flagD = widthD*heightD;
viscousDampingRatio = (1/(4*pi)) * (flagD/linearD);

widthD2 = x3-x2;
heightD2 = y2-y1;
flagD2 = widthD2*heightD2;
viscousDampingRatio2 = (1/(4*pi)) * (flagD2/linearD);

%% Loop!

g = 9.81;
R = 4;
Tn = 0.7;
wn = 2*pi/Tn;
Kinitial = wn*wn*MAll*0.8;
WAll = MAll*g;
myAlpha = 0.8;
myGamma = 0.1;
NN = 0;
NPTpers = length(PTpers);
NEQs = 60;
VDRmatrix = zeros(NPTpers,NEQs);

for ii = 1:NPTpers
    
    for jj = 1:NEQs
        
        NN = NN + 1;
        
        Gup = Gups(jj);
        PTper = PTpers(ii);
        Mover = Gup*g*MAll*HM/R;
        MED = Mover*(100-PTper)/100;
        MPT = Mover - MW - MED;
        initialFPT = MPT/B;
        maxD = dispx(NN,6);
        maxF = max([reactx(NN,5) reactx(NN,6)]);
        linearD = maxF*maxD/2;
        
        y1 = (MW + MPT)/H;
        x1 = y1/Kinitial;
        
        y2 = Mover/H;
        x2 = x1 + (MED/H)/(Kinitial*myAlpha);
        
        y3 = maxF;
        x3 = maxD;
        
        widthD = x3-x2;
        heightD = y2-y1;
        flagD = widthD2*heightD2;
        viscousDampingRatio = (1/(4*pi)) * (flagD2/linearD);
        
        VDRmatrix(ii,jj) = viscousDampingRatio;
        
    end
    
end


%% Height / period and deflections from equivalent static method ??
close all;

siteSpectrum = load('siteSpectrum.txt');

T1 = testDetails(2);
R = testDetails(8);
% AR = 4;

Mt = 21e5;
g = 9.81;
Wt = Mt*g;
NT1 = int8(T1/0.05);

K1 = (2*pi/T1)^2*Mt;
CT1 = siteSpectrum(NT1);
mu = 4;
Sp = 0.7;
dispFactor = 3/2;

if T1 < 0.7
    Kmu = (mu-1)*T1/0.7 + 1;
else
    Kmu = mu;
end

CdT1 = CT1*Sp./Kmu;
V = CdT1*Wt;
deltaY = V/K1*dispFactor;
delta = 1.3*mu*deltaY;
deltaMod = 3/2*delta;

% B = 2.5;
% R = 4;
MED = 0.34;
myAlpha = 0.8;
myBeta = 0.1;
% H = AR*2*B; % H = T1^(1/0.75)./0.0625;
 % H/10;

VED = CT1*Wt/R;
Vuplift = (1-MED)*VED;
Kuplift = myAlpha*K1;
KED = myBeta*K1;
deltaUplift = Vuplift/K1*dispFactor;
deltaED = deltaUplift + (VED-Vuplift)/Kuplift*dispFactor;
Vmax = VED + KED*(delta-deltaED);
deltaElastic = Vmax/K1;
deltaRotation = delta-deltaElastic;
theta = deltaRotation/H;
xuplift = 2*B.*theta;

% dispArray = [0 deltaUplift deltaED delta deltaRotation];
% FArray = [0 Vuplift VED Vmax 0];

% Complete the flag
VflagSE = Vmax-(VED-Vuplift);
dispflagSE = delta - (Vmax-VflagSE)/K1;

dispArray = [0 deltaUplift deltaED delta dispflagSE deltaUplift];
FArray = [0 Vuplift VED Vmax VflagSE Vuplift];

% plot result

% figure
% plot(dispArray,FArray,'o--')
% grid on
% xlabel('Displacement')
% ylabel('Force')

figure
plot(dispx(:,10)-dispx(:,1),[-reactx(:,5) -reactx(:,6)])
hold on
plot(dispArray,FArray,'go--','lineWidth',2)
grid on
legend('LHS rocking edge','RHS rocking edge','ESM','Location','Best')
xlabel('Displacement (m)')
ylabel('Force at rocking edge (N)')
title('Deflection at top of frame')
fileLoc = strcat(filedir,'Figures/EQ',num2str(EQnum),'T',num2str(T1*10),'AR',num2str(AR));
saveas(gcf,fileLoc,'meta')
saveas(gcf,fileLoc,'fig')

return

%% Sample results figures

% MATERIAL MODELS

% Horizontal support - left
figure
plot(dispx(:,3)*1000,-reactx(:,5)/1e6,'b')
grid on
xlabel('Displacement (mm)')
ylabel({'Node 5 Horizontal','Force (MN)'})
axis([-inf inf -inf 1])
xt = get(gca, 'XTick');
set(gca, 'FontSize', 34)
set(gca, 'FontWeight', 'bold')
fileLoc = 'Figures/SupportLHSlargeprint';
saveas(gcf,fileLoc,'meta')
saveas(gcf,fileLoc,'fig')

figure
plot(dispx(:,3)*1000,-reactx(:,5)/1e6,'b')
grid on
xlabel('Displacement (mm)')
ylabel('Node 5 Horizontal Force (MN)')
axis([-inf inf -inf 1])
xt = get(gca, 'XTick');
set(gca, 'FontSize', 22)
fileLoc = 'Figures/SupportLHS';
saveas(gcf,fileLoc,'meta')
saveas(gcf,fileLoc,'fig')

% Horizontal support - right
figure
plot(dispx(:,4)*1000,-reactx(:,6)/1e6,'r')
grid on
xlabel('Displacement (mm)')
ylabel({'Node 6 Horizontal','Force (MN)'})
axis([-inf inf -1 inf])
xt = get(gca, 'XTick');
set(gca, 'FontSize', 34)
set(gca, 'FontWeight', 'bold')
fileLoc = 'Figures/SupportRHSlargeprint';
saveas(gcf,fileLoc,'meta')
saveas(gcf,fileLoc,'fig')

figure
plot(dispx(:,4)*1000,-reactx(:,6)/1e6,'r')
grid on
xlabel('Displacement (mm)')
ylabel('Node 6 Horizontal Force (MN)')
axis([-inf inf -1 inf])
xt = get(gca, 'XTick');
set(gca, 'FontSize', 22)
fileLoc = 'Figures/SupportRHS';
saveas(gcf,fileLoc,'meta')
saveas(gcf,fileLoc,'fig')

% Device and rocking edge - left
figure
plot(dispy(:,3)*1000,-reacty(:,5)/1e6,'b')
grid on
axis([-5 55 -30 5])
xlabel('Displacement (mm)')
ylabel({'Node 5 Vertical','Force (MN)'})
xt = get(gca, 'XTick');
set(gca, 'FontSize', 34)
set(gca, 'FontWeight', 'bold')
fileLoc = 'Figures/GNGLHSlargeprint';
saveas(gcf,fileLoc,'meta')
saveas(gcf,fileLoc,'fig')

figure
plot(dispy(:,3)*1000,-reacty(:,5)/1e6,'b')
grid on
axis([-5 55 -30 5])
xt = get(gca, 'XTick');
set(gca, 'FontSize', 22)
xlabel('Displacement (mm)')
ylabel('Node 5 Vertical Force (MN)')
fileLoc = 'Figures/GNGLHS';
saveas(gcf,fileLoc,'meta')
saveas(gcf,fileLoc,'fig')

figure
plot(dispy(:,3)*1000,-reacty(:,5),'b')
grid on
xlabel('Displacement (mm)')
ylabel('Node 5 Vertical Force (MN)')
axis([0 55 0 5e6])
xt = get(gca, 'XTick');
set(gca, 'FontSize', 22)
fileLoc = 'Figures/GNGLHScloseup';
saveas(gcf,fileLoc,'meta')
saveas(gcf,fileLoc,'fig')

% Device and rocking edge - right
figure
plot(dispy(:,4)*1000,-reacty(:,6)/1e6,'r')
grid on
axis([-5 50 -30 5])
xlabel('Displacement (mm)')
ylabel({'Node 6 Vertical','Force (MN)'})
xt = get(gca, 'XTick');
set(gca, 'FontSize', 34)
set(gca, 'FontWeight', 'bold')
fileLoc = 'Figures/GNGRHSlargeprint';
saveas(gcf,fileLoc,'meta')
saveas(gcf,fileLoc,'fig')

% Post-tensioning
figure
plot(dispy(:,1)*1000,-reacty(:,7)/1e6)
grid on
axis([-1 30 12 30])
xlabel('Displacement (mm)')
ylabel({'Node 7 Vertical','Force (MN)'})
xt = get(gca, 'XTick');
set(gca, 'FontSize', 34)
set(gca, 'FontWeight', 'bold')
fileLoc = 'Figures/PTlargeprint';
saveas(gcf,fileLoc,'meta')
saveas(gcf,fileLoc,'fig')

figure
plot(dispy(:,1)*1000,-reacty(:,7)/1e6)
grid on
axis([-1 30 12 30])
xt = get(gca, 'XTick');
set(gca, 'FontSize', 22)
xlabel('Displacement (mm)')
ylabel('Node 7 Vertical Force (MN)')
fileLoc = 'Figures/PT';
saveas(gcf,fileLoc,'meta')
saveas(gcf,fileLoc,'fig')

% RESPONSE

% Drift at top of frame
figure
plot(time,(dispx(:,10)-dispx(:,1))/H*100)
grid on
xlabel('Time (s)')
ylabel('Drift at node 10 (%)')
xt = get(gca, 'XTick');
set(gca, 'FontSize', 22)
fileLoc = 'Figures/Drift';
saveas(gcf,fileLoc,'meta')
saveas(gcf,fileLoc,'fig')

% Uplift at rocking edges
figure
plot(time,[dispy(:,3),dispy(:,4)]*1000)
grid on
ylim([-1 60])
xlabel('Time (s)')
ylabel('Uplift (mm)')
legend('LHS rocking edge','RHS rocking edge')
xt = get(gca, 'XTick');
set(gca, 'FontSize', 22)
fileLoc = 'Figures/Uplift';
saveas(gcf,fileLoc,'meta')
saveas(gcf,fileLoc,'fig')

% Hysteresis
figure
plot((dispx(:,10)-dispx(:,1))/H*100,[-reactx(:,5) -reactx(:,6)]/1e6)
grid on
axis([-1.5 1.5 -8 8])
xlabel('Drift at Node 10 (%)')
ylabel('Rocking base node horizontal force (MN)')
legend('LHS rocking edge','RHS rocking edge','Location','Best')
xt = get(gca, 'XTick');
set(gca, 'FontSize', 22)
fileLoc = 'Figures/FD';
saveas(gcf,fileLoc,'meta')
saveas(gcf,fileLoc,'fig')

% Dissipater plastic demand
figure
plot(time,[demand(:,1) demand(:,2)]*1000)
grid on
xlabel('Time (s)')
ylabel('Inelastic displacement (mm)')
legend('LHS device','RHS device','Location','Best')
xt = get(gca, 'XTick');
set(gca, 'FontSize', 22)
fileLoc = 'Figures/Demand';
saveas(gcf,fileLoc,'meta')
saveas(gcf,fileLoc,'fig')

% Device ratchet count
figure
plot(time,[ratchet(:,1) ratchet(:,2)])
grid on
xlabel('Time (s)')
ylabel('Number of ratcheting actions')
legend('LHS device','RHS device','Location','Best')
xt = get(gca, 'XTick');
set(gca, 'FontSize', 22)
fileLoc = 'Figures/RatchetCount';
saveas(gcf,fileLoc,'meta')
saveas(gcf,fileLoc,'fig')

% GROUND MOTION

scaleFactors = load('scaleFactors.txt');
AllEQs = load('earthquake signals/Standard Set/t_0_001s/EQ_records_All.mat');
timeGM = [0 AllEQs.Time_vector];
Accels = AllEQs.Accel_matrix;
NGM = size(Accels,2);
Accels = [zeros(1,NGM); Accels];
Ts = (0.05:0.05:5);
NT = length(Ts);
GMno = 1;
thePeriod = testDetails(2);
TnRef = find(Ts == thePeriod);
TnRef = 14;

figure
plot(timeGM,scaleFactors(TnRef,GMno)*Accels(:,GMno))
grid on
% axis([-inf inf -3 3])
xlabel('Time (s)')
ylabel('Acceleration (ms^-^2)')
xt = get(gca, 'XTick');
set(gca, 'FontSize', 34)
set(gca, 'FontWeight', 'bold')
fileLoc = 'Figures/GM01T07largeprint';
saveas(gcf,fileLoc,'meta')
saveas(gcf,fileLoc,'fig')

figure
plot(timeGM,scaleFactors(TnRef,GMno)*Accels(:,GMno))
grid on
xlabel('Time (s)')
ylabel('Acceleration (ms^-^2)')
xt = get(gca, 'XTick');
set(gca, 'FontSize', 22)
fileLoc = 'Figures/GM01T07';
saveas(gcf,fileLoc,'meta')
saveas(gcf,fileLoc,'fig')

% FANCY COMMON TIME AXIS PLOT

figure('Position',[50 200 1200 800])
ax1 = subplot(3,1,1);
plot(timeGM,scaleFactors(TnRef,GMno)*Accels(:,GMno))
grid on
ylim([-4 4])
xlim([0 50])
ylabel({'Acceleration','(ms^-^2)'})
xt = get(gca, 'XTick');
set(gca, 'FontSize', 22)
ax2 = subplot(3,1,2);
plot(time,(dispx(:,10)-dispx(:,1))/H*100)
grid on
ylim([-1.5 1.5])
xlim([0 50])
ylabel({'Drift at node 10','(%)'})
xt = get(gca, 'XTick');
set(gca, 'FontSize', 22)
ax3 = subplot(3,1,3);
plot(time,[dispy(:,3),dispy(:,4)]*1000)
grid on
ylim([-inf 60])
xlim([0 50])
ylabel({'Uplift','(mm)'})
xlabel('Time (s)')
legend('LHS rocking edge','RHS rocking edge')
xt = get(gca, 'XTick');
set(gca, 'FontSize', 22)
linkaxes([ax1,ax2,ax3],'x')
fileLoc = 'Figures/GM01DriftUplift';
saveas(gcf,fileLoc,'meta')
saveas(gcf,fileLoc,'fig')


%% logarithmic decrement stuff
figure(101), figure(102), figure(103)
close([101 102 103])

% figure
% plot(time,[dispy(:,3),dispy(:,4)],'linewidth',1.5)
% grid on
% xlabel('Time')
% ylabel('Displacement')
% title('RE vertical disp.')
% 
% figure
% plot(time,dispx(:,10))
% grid on
% xlabel('Time')
% ylabel('Displacement')
% title('Roof displacement')

% get the damping ratio from the logarithmic decrement

x = dispx(:,10);
% x = [dispy(:,3) -dispy(:,4)];
% t = time;
t = [time time];

% find peaks
diff_sign_diff = (diff(sign(diff(x))));
neg_peaks = find(diff_sign_diff == 2)+1;
pos_peaks = find(diff_sign_diff == -2)+1;

[pks,neg_peaks] = findpeaks(-x);
[pks,pos_peaks] = findpeaks(x);

figure(101)
plot(t,x);grid on;hold on
plot(t(pos_peaks),x(pos_peaks),'ro');
plot(t(neg_peaks),x(neg_peaks),'go');

% check the damping ratio from the logarithmic decrement using the linear
% simulation
pos_ratios = x(pos_peaks(2:end))./x(pos_peaks(1:end-1));
neg_ratios = x(neg_peaks(2:end))./x(neg_peaks(1:end-1));

% my_pos_peaks = [0.1234 0.07708 0.05507 0.04112 0.03135 0.02412 0.01862 0.01434 0.01096];
% my_neg_peaks = [0.09449 0.06463 0.04742 0.03584 0.02748 0.0212 0.01634 0.01255];
% 
% pos_ratios = my_pos_peaks(2:end)./my_pos_peaks(1:end-1);
% neg_ratios = my_neg_peaks(2:end)./my_neg_peaks(1:end-1);

figure(102)
plot(pos_ratios,'bo-'),hold on
plot(neg_ratios,'ro-')
grid on

overall_mean_ratio = (mean(pos_ratios)+mean(neg_ratios))/2
plot(overall_mean_ratio*ones(size(pos_ratios)),'m')

identified_damping_ratio = 1/sqrt(1+(2*pi/log(1/overall_mean_ratio))^2)

peakRats = 0.01:0.01:1;
logDecs = log(1./peakRats);
dampRats = 1./sqrt(1 + (2*pi./logDecs).^2);

figure(103)
plot(peakRats,dampRats)
grid on
xlabel('Peak ratio')
ylabel('Damping ratio')

%%

% figure
% plot(time(2:end),[KsLHS KsRHS])
% grid on
% xlabel('Time (s)')
% ylabel('Stiffness (N/m)') 
% axis([-inf inf -2*Ktrue 2*Ktrue])


%% Acceleration ramp figures

% theEnd = size(dispx,1);
% 
% figure
% plot(AccelRamp(1:theEnd),dispx(:,2))
% grid on
% xlabel('Acceleration (m/s/s)')
% ylabel('Displacement (m)')
% title('Displacement at mass')
% 
% figure
% plot(AccelRamp(1:theEnd),reactx(:,5))
% grid on
% xlabel('Acceleration (m/s/s)')
% ylabel('Base force (N)')



figure
plot(time,[demand(:,1) demand(:,2)])
grid on
xlabel('Time, s')
ylabel('Inelastic displacement, m')
title('Cumulative inelastic dissipater demand')

figure
plot(time,[ratchet(:,1) ratchet(:,2)])
grid on
xlabel('Time, s')
ylabel('No. ratcheting actions')
title('GnG ratchet count')


figure
plot(dispx(:,3),dispy(:,3))
grid on
axis([-inf inf min(dispx(:,3)) max(dispx(:,3))])
xlabel('X displacement')
ylabel('Y displacement')
title('Left rocking edge motion')

figure
plot(dispx(:,4),dispy(:,4))
grid on
axis([-inf inf min(dispx(:,4)) max(dispx(:,4))])
xlabel('X displacement')
ylabel('Y displacement')
title('Right rocking edge motion')


figure
plot(time,[-reactx(:,5)-reactx(:,6)])
grid on
figure
plot(time,dispx(:,2))
grid on

%%
%------------dissipater elongation with random load reversals------------%
testTime = 0:0.005:5;
testData = 0.03*sin(testTime) + 0.005*sin(12*testTime);
figure
plot(testData)
grid on
testyield = 8e-3;
pks = findpeaks(testData);
pksloc = pks > testyield;
bigpks = pks(pksloc);
elong = bigpks - testyield;
totalelong = sum(elong);
fprintf('Test elongation = %5.2f units\n',totalelong*1000)

N = length(testData);
elongtest = 0;
for i = 2:N
    % LHS loop
    if testData(i) > testyield
        if testData(i) > testData(i-1)
            elongtest = elongtest + (testData(i) - testData(i-1));
        end
    end
end
fprintf('Test elongation = %5.2f units\n',elongtest*1000)

% [invtrs,trslocs] = findpeaks(-testData);
% trs = testData(trslocs);
% bigtrs = trs > testyield;
% 
% invertedtestData = max(testData(:)) - testData;
% trs = findpeaks(invertedtestData);
% 
% trsloc = islocalmin(testData);
% trs = testData(trsloc);



N = length(disp3y);
elongLHS = 0;
elongRHS = 0;
for i = 2:N
    % LHS loop
    if dispy(i,3) > dyield
        if dispy(i,3) > dispy(i-1,3)
            elongLHS = elongLHS + (dispy(i,3) - dispy(i-1,3));
        end
    end
    % RHS loop
    if dispy(i,4) > dyield
        if dispy(i,4) > dispy(i-1,4)
            elongRHS = elongRHS + (dispy(i,4) - dispy(i-1,4));
        end
    end
end
fprintf('Elongation of LHS dissipater = %5.2f mm\n',elongLHS*1000)
fprintf('Elongation of RHS dissipater = %5.2f mm\n',elongRHS*1000)

%------------dissipater elongation with random load reversals------------%

return

% Rocking edge forces - can clearly see weight force evenly distributed
% initially and then zero force in one edge while full weight force in the
% other edge during uplift and crossover during transitions
figure
plot(time,[reacty(:,5),reacty(:,6)])
grid on
legend('left foot','right foot')
xlabel('time')
ylabel('vertical reaction force at rocking edge')

figure
plot(time,dispx(:,2))
grid on

%% Animation of results
close all
% aspect ratio
H = 5; % total height
HM = H*2/3; % mass height
b = 2.5; % half_width
L = 3*b; % leaning column offset
N = length(time);
filename = 'rockingVid.gif';
saveitpls = 0; % set to 1 to save gif

% m = 1e7;
% g = 9.81;
% k = 200e9;
% k2 = 200e9*0.1/H;
% 
% uc = m*g*b/k/H;
% uc2 = m*g*b/k2/H;

% undeflected node positions
x01 = 0;
y01 = 0;
x02 = 0;
y02 = HM;
x03 = -b;
y03 = 0;
x04 = b;
y04 = 0;
x09 = L;
y09 = HM;
x010 = 0;
y010 = H;

% deflected node positions
pos1x = x01 + dispx(:,1);
pos1y = y01 + dispy(:,1);
pos2x = x02 + dispx(:,2);
pos2y = y02 + dispy(:,2);
pos3x = x03 + dispx(:,3);
pos3y = y03 + dispy(:,3);
pos4x = x04 + dispx(:,4);
pos4y = y04 + dispy(:,4);
pos9x = x09 + dispx(:,9);
pos9y = y09 + dispy(:,9);
pos10x = x010 + dispx(:,10);
pos10y = y010 + dispy(:,10);
% top undeflected
hyp = 2*b;
opp = pos3y - pos4y; 
theta = asin(opp/hyp);
pos10xun = pos1x + H*sin(theta);
pos10yun = pos1y + H*cos(theta);

myfig = figure('Position',[100 100 850 800]);
hline11 = plot(NaN,'ko-','linewidth',3,'MarkerSize',12);
hold on
hline12 = plot(NaN,'ro--','linewidth',3,'MarkerSize',12);
hline13 = plot(NaN,'go-','linewidth',3,'MarkerSize',12);
legend('Flexible frame','Reference','Leaning column','Location','Best')
% axis([-2*b L+b -0.1*H 1.1*H])
axis([-2*b L+b -2*b L+b])
grid on

myfig2 = figure('Position',[1000 100 850 300]);
hline21 = plot(NaN,'ko-','linewidth',3,'MarkerSize',12);
hold on
hline22 = plot(NaN,'ro--','linewidth',3,'MarkerSize',12);
legend('Flexible frame','Reference','Location','SouthEast')
madMax = max([max(dispy(:,3)) max(dispy(:,4))]);
axis([-1.1*b 1.1*b -madMax/2 madMax])
grid on

xloc22 = [x03,x01,x04];
yloc22 = [y03,y01,y04];
set(hline22,'XData',xloc22)
set(hline22,'YData',yloc22)

myfig3 = figure('Position',[1000 500 850 400]);
subplot(1,2,1)
hline31 = plot(NaN,'b-','linewidth',3,'MarkerSize',12);
hold on
hline32 = plot(NaN,'r-','linewidth',3,'MarkerSize',12);
axis([0 time(end) 0 madMax])
grid on
xlabel('Time')
ylabel('Displacement')
title('Rocking edge displacement')
legend('Left','Right','Location','NorthEast')

subplot(1,2,2)
hline41 = plot(NaN,'g-','linewidth',3,'MarkerSize',12);
axis([0 time(end) min(dispx(:,2)) max(dispx(:,2))])
grid on
xlabel('Time')
ylabel('Displacement')
title('Roof displacement')

for kk = 1:100:N
    
    % Figure 1 - system
    % deflected positions
    xloc11 = [pos3x(kk),pos4x(kk),pos1x(kk),pos2x(kk),pos10x(kk)];
    yloc11 = [pos3y(kk),pos4y(kk),pos1y(kk),pos2y(kk),pos10y(kk)];
    set(hline11,'XData',xloc11)
    set(hline11,'YData',yloc11)
    % undeflected vertical line
    xloc12 = [pos1x(kk),pos10xun(kk)];
    yloc12 = [pos1y(kk),pos10yun(kk)];
    set(hline12,'XData',xloc12)
    set(hline12,'YData',yloc12)
    % leaning column
    xloc13 = [L,pos9x(kk),pos2x(kk)];
    yloc13 = [0,pos9y(kk),pos2y(kk)];
    set(hline13,'XData',xloc13)
    set(hline13,'YData',yloc13)
    
    % Figure 2 - base
    % base
    xloc21 = [pos3x(kk),pos1x(kk),pos4x(kk)];
    yloc21 = [pos3y(kk),pos1y(kk),pos4y(kk)];
    set(hline21,'XData',xloc21)
    set(hline21,'YData',yloc21)
    
    % Figure 3 - displacements
    % rocking edges
    x31 = time(1:kk);
    y31 = dispy(1:kk,3);
    set(hline31,'XData',x31)
    set(hline31,'YData',y31)
    x32 = time(1:kk);
    y32 = dispy(1:kk,4);
    set(hline32,'XData',x32)
    set(hline32,'YData',y32)
    % roof
    x41 = time(1:kk);
    y41 = dispx(1:kk,2);
    set(hline41,'XData',x41)
    set(hline41,'YData',y41)
    
    drawnow
    
    pause(0.1)
    
    fprintf('Time = %d\n',time(kk))
    
    % Capture the plot as an image
    frame = getframe(myfig);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    % Write to the GIF File
    if saveitpls == 1
        if kk == 1
            imwrite(imind,cm,filename,'gif', 'Loopcount',1);
        else
            imwrite(imind,cm,filename,'gif','WriteMode','append');
        end
    end
    
end

%% Animation of results: force - displacement
% close all

N = length(time);
myfig = figure('Position',[100 200 600 500]);
hline11 = plot(NaN,'b-');
hold on
hline12 = plot(NaN,'r-');
hline13 = plot(NaN,'go','lineWidth',2,'MarkerSize',8);
hline14 = plot(NaN,'go','lineWidth',2,'MarkerSize',8);
legend('Left foot','Right foot','Location','Best')
axis([min(dispx(:,2)) max(dispx(:,2)) min(-reactx(:,5)) max(-reactx(:,6))])
grid on
xlabel('Lateral deflection')
ylabel('Base reaction force')
title('System force-displacement')

myfig2 = figure('Position',[700 200 600 500]);
hline21 = plot(NaN,'k-');
hold on
hline22 = plot(NaN,'go','lineWidth',2,'MarkerSize',8);
axis([min(dispy(:,3)) max(dispy(:,3)) min(-reacty(:,5)) max(-reacty(:,5))])
grid on
xlabel('Displacement')
ylabel('Reaction force')
title('LHS GNG - vertical')

myfig3 = figure('Position',[1300 200 600 500]);
hline31 = plot(NaN,'k-');
hold on
hline32 = plot(NaN,'go','lineWidth',2,'MarkerSize',8);
axis([min(dispy(:,1)) max(dispy(:,1)) min(-reacty(:,7)) max(-reacty(:,7))])
grid on
xlabel('Displacement')
ylabel('Reaction force')
title('post-tensioning')


for kk = 4200:1:N
    
    % System force - displacement
    xleft = dispx(1:kk,2);
    yleft = -reactx(1:kk,5);
    set(hline11,'XData',xleft)
    set(hline11,'YData',yleft)
    xright = dispx(1:kk,2);
    yright = -reactx(1:kk,6);
    set(hline12,'XData',xright)
    set(hline12,'YData',yright)
    set(hline13,'XData',dispx(kk,2))
    set(hline13,'YData',-reactx(kk,5))
    set(hline14,'XData',dispx(kk,2))
    set(hline14,'YData',-reactx(kk,6))
    
    set(hline21,'XData',dispy(1:kk,3))
    set(hline21,'YData',-reacty(1:kk,5))
    set(hline22,'XData',dispy(kk,3))
    set(hline22,'YData',-reacty(kk,5))
    
    set(hline31,'XData',dispy(1:kk,1))
    set(hline31,'YData',-reacty(1:kk,7))
    set(hline32,'XData',dispy(kk,1))
    set(hline32,'YData',-reacty(kk,7))
    
    drawnow
    
%     pause(0.05)
    
    if (rem(time(kk),0.1) < 1e-4) 
        fprintf('Time = %d\n',time(kk))
    end
    
end
