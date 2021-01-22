%% Load OpenSees data from sample parameter study with 5 screens

% clear previous work
close all; clear; clc;

%% Load in the data

% geometry
B = 2.5;
PTper = 66;

% file directory
filedir = 'finalResults3_10_11_R4/';

% load the data
theData1 = load(strcat(filedir,'testResults.txt'));

% Input arrays
gMotionsRaw = theData1(:,1);
TnsRaw      = theData1(:,2);
SFsRaw      = theData1(:,3);
ARsRaw      = theData1(:,4);
PitchsRaw   = theData1(:,5);
PTpersRaw   = theData1(:,6);
Nanalyses = numel(TnsRaw);

% Output arrays
Fyields     = theData1(:,7);
oks         = theData1(:,8);
eigs        = theData1(:,9);
demandLHS   = theData1(:,10);
demandRHS   = theData1(:,11);
ratsLHS     = theData1(:,12);
ratsRHS     = theData1(:,13);

% Lists of unique input values
gMotions = unique(gMotionsRaw);
Tns = unique(TnsRaw);
ARs = unique(ARsRaw);
Pitchs = unique(PitchsRaw);
PTpers = unique(PTpersRaw);

TnsAR2 = [0.3 0.4];
TnsAR4 = [0.5 0.6 0.7 0.8];
TnsAR6 = [0.8 0.9 1.0 1.1 1.2];
TnsAR8 = [1.0 1.1 1.2 1.3 1.4 1.5 1.6];
TnsCell = {TnsAR2 TnsAR4 TnsAR6 TnsAR8};

gMotions3 = unique(gMotionsRaw,'stable');

wns = 2*pi./Tns;

% lengths of unique input arrays
NgMotions = length(gMotions);
NTns = length(Tns);
NARs = length(ARs);
NPitchs = length(Pitchs);
NPTpers = length(PTpers);
NTnsArray = [2 4 5 7];
NTns2 = max(NTnsArray);

% Node recorders
DispData        = load(strcat(filedir,'nodeD.txt'));
AccelData       = load(strcat(filedir,'nodeA.txt'));
ReactionData    = load(strcat(filedir,'nodeR.txt'));

dispx   = DispData(:,1:3:end);
dispy   = DispData(:,2:3:end);
dispz   = DispData(:,3:3:end);

Ndispx = size(dispx,2);

% Accelerations
accelx   = AccelData(:,1:2:end);
accely   = AccelData(:,2:2:end);
Naccelx = size(accelx,2);

% reactions
reactx  = ReactionData(:,1:2:end);
reacty  = ReactionData(:,2:2:end);
Nreactx = size(reactx,2);

%% Preallocation

% Preallocate matrices (not for geometric means, NaNs for unstable results)

SAMatrixMass = zeros(NgMotions,NARs,NTns2,NPitchs);
SAMatrixTop = zeros(NgMotions,NARs,NTns2,NPitchs);
HMatrix = zeros(NARs,NTns2,NPitchs);
HMatrixAll = zeros(NgMotions,NARs,NTns2,NPitchs);
PitchMatrix = zeros(NgMotions,NARs,NTns2,NPitchs);
FyieldMatrix = zeros(NARs,NTns2,NPitchs);
okMatrix = zeros(NgMotions,NARs,NTns2,NPitchs);
eigsMatrix = zeros(NgMotions,NARs,NTns2,NPitchs);
demandLHSMatrix = zeros(NgMotions,NARs,NTns2,NPitchs);
demandRHSMatrix = zeros(NgMotions,NARs,NTns2,NPitchs);
ratsLHSMatrix = zeros(NgMotions,NARs,NTns2,NPitchs);
ratsRHSMatrix = zeros(NgMotions,NARs,NTns2,NPitchs);

dempersLHSMatrix = ones(NgMotions,NARs,NTns2,NPitchs);
dempersRHSMatrix = ones(NgMotions,NARs,NTns2,NPitchs);

dispxMatrix = zeros(NgMotions,Ndispx,NARs,NTns2,NPitchs);
dispyMatrix = zeros(NgMotions,Ndispx,NARs,NTns2,NPitchs);
dispzMatrix = zeros(NgMotions,Ndispx,NARs,NTns2,NPitchs);
accelxMatrix = zeros(NgMotions,Naccelx,NARs,NTns2,NPitchs);
accelyMatrix = zeros(NgMotions,Naccelx,NARs,NTns2,NPitchs);
reactxMatrix = zeros(NgMotions,Nreactx,NARs,NTns2,NPitchs);
reactyMatrix = zeros(NgMotions,Nreactx,NARs,NTns2,NPitchs);

% Preallocate matrices for geometric means (no NaNs, using ones to allow
% for nth root calculation, '2' tag)

demandLHSMatrix2 = ones(NgMotions,NARs,NTns2,NPitchs);
demandRHSMatrix2 = ones(NgMotions,NARs,NTns2,NPitchs);
ratsLHSMatrix2 = ones(NgMotions,NARs,NTns2,NPitchs);
ratsRHSMatrix2 = ones(NgMotions,NARs,NTns2,NPitchs);

dempersLHSMatrix2 = ones(NgMotions,NARs,NTns2,NPitchs);
dempersRHSMatrix2 = ones(NgMotions,NARs,NTns2,NPitchs);

dispxMatrix2 = ones(NgMotions,Ndispx,NARs,NTns2,NPitchs);
dispyMatrix2 = ones(NgMotions,Ndispx,NARs,NTns2,NPitchs);
dispzMatrix2 = ones(NgMotions,Ndispx,NARs,NTns2,NPitchs);
accelxMatrix2 = ones(NgMotions,Naccelx,NARs,NTns2,NPitchs);
accelyMatrix2 = ones(NgMotions,Naccelx,NARs,NTns2,NPitchs);
reactxMatrix2 = ones(NgMotions,Nreactx,NARs,NTns2,NPitchs);
reactyMatrix2 = ones(NgMotions,Nreactx,NARs,NTns2,NPitchs);

% Number of stable results
NMatrix = zeros(NARs,NTns2,NPitchs);

%% Create matrices

ctr = 1;

for aa = 1:NgMotions
    
    GMref = gMotions3(aa);
    
    for ii = 1:NARs
        
        for jj = 1:NTnsArray(ii)
            
            for kk = 1:NPitchs
                
                okMatrix(GMref,ii,jj,kk) = oks(ctr);
                eigsMatrix(GMref,ii,jj,kk) = eigs(ctr);
                PitchMatrix(aa,ii,jj,kk) = Pitchs(kk);
                
                thisH = ARs(ii)*2*B;
                HMatrix(ii,jj,kk) = thisH;
                HMatrixAll(GMref,ii,jj,kk) = thisH;
                
                
                if oks(ctr) ~= 0
                    
                    SAMatrixMass(GMref,ii,jj,kk) = NaN;
                    SAMatrixTop(GMref,ii,jj,kk) = NaN;
                    
                    FyieldMatrix(GMref,ii,jj,kk) = NaN;
                    demandLHSMatrix(GMref,ii,jj,kk) = NaN;
                    demandRHSMatrix(GMref,ii,jj,kk) = NaN;
                    ratsLHSMatrix(GMref,ii,jj,kk) = NaN;
                    ratsRHSMatrix(GMref,ii,jj,kk) = NaN;
                    
                    dispxMatrix(GMref,:,ii,jj,kk) = NaN;
                    dispyMatrix(GMref,:,ii,jj,kk) = NaN;
                    dispzMatrix(GMref,:,ii,jj,kk) = NaN;
                    accelxMatrix(GMref,:,ii,jj,kk) = NaN;
                    accelyMatrix(GMref,:,ii,jj,kk) = NaN;
                    reactxMatrix(GMref,:,ii,jj,kk) = NaN;
                    reactyMatrix(GMref,:,ii,jj,kk) = NaN;
                    
                    dempersLHSMatrix(GMref,ii,jj,kk) = NaN;
                    dempersRHSMatrix(GMref,ii,jj,kk) = NaN;
                    
                else
                    
                    SAMatrixMass(GMref,ii,jj,kk) = (wns(ii)^2)*dispx(ctr,2);
                    SAMatrixTop(GMref,ii,jj,kk) = (wns(ii)^2)*dispx(ctr,6);
                    
                    NMatrix(ii,jj,kk) = NMatrix(ii,jj,kk) + 1;
                    
                    FyieldMatrix(GMref,ii,jj,kk) = Fyields(ctr);
                    demandLHSMatrix(GMref,ii,jj,kk) = demandLHS(ctr);
                    demandRHSMatrix(GMref,ii,jj,kk) = demandRHS(ctr);
                    ratsLHSMatrix(GMref,ii,jj,kk) = ratsLHS(ctr);
                    ratsRHSMatrix(GMref,ii,jj,kk) = ratsRHS(ctr);
                    
                    dispxMatrix(GMref,:,ii,jj,kk) = dispx(ctr,:);
                    dispyMatrix(GMref,:,ii,jj,kk) = dispy(ctr,:);
                    dispzMatrix(GMref,:,ii,jj,kk) = dispz(ctr,:);
                    accelxMatrix(GMref,:,ii,jj,kk) = accelx(ctr,:);
                    accelyMatrix(GMref,:,ii,jj,kk) = accely(ctr,:);
                    reactxMatrix(GMref,:,ii,jj,kk) = reactx(ctr,:);
                    reactyMatrix(GMref,:,ii,jj,kk) = reacty(ctr,:);
                    
                    demandLHSMatrix2(GMref,ii,jj,kk) = demandLHS(ctr);
                    demandRHSMatrix2(GMref,ii,jj,kk) = demandRHS(ctr);
                    ratsLHSMatrix2(GMref,ii,jj,kk) = ratsLHS(ctr);
                    ratsRHSMatrix2(GMref,ii,jj,kk) = ratsRHS(ctr);
                    
                    dispxMatrix2(GMref,:,ii,jj,kk) = dispx(ctr,:);
                    dispyMatrix2(GMref,:,ii,jj,kk) = dispy(ctr,:);
                    dispzMatrix2(GMref,:,ii,jj,kk) = dispz(ctr,:);
                    accelxMatrix2(GMref,:,ii,jj,kk) = accelx(ctr,:);
                    accelyMatrix2(GMref,:,ii,jj,kk) = accely(ctr,:);
                    reactxMatrix2(GMref,:,ii,jj,kk) = reactx(ctr,:);
                    reactyMatrix2(GMref,:,ii,jj,kk) = reacty(ctr,:);
                    
                    if kk ~= 1
                        
                        if (demandLHS(ctr)*demandLHSMatrix(GMref,ii,jj,1)) ~= 0
                            
                            dempersLHSMatrix2(GMref,ii,jj,kk) = demandLHS(ctr)...
                                /demandLHSMatrix(GMref,ii,jj,1);
                            
                        end
                        
                        if (demandRHS(ctr)*demandRHSMatrix(GMref,ii,jj,1)) ~= 0
                            
                            dempersRHSMatrix2(GMref,ii,jj,kk) = demandRHS(ctr)...
                                /demandRHSMatrix(GMref,ii,jj,1);
                            
                        end
                        
                    end
                    
                end
                
                ctr = ctr + 1;
                
            end
            
        end
    end
    
end

dempersLHSMatrix(:,:,:,1) = 1;
dempersRHSMatrix(:,:,:,1) = 1;

dempersLHSMatrix2(:,:,:,1) = 1;
dempersRHSMatrix2(:,:,:,1) = 1;

%% DIGGING INTO STABILITY

Nfails = sum(okMatrix(:))/-3;
fprintf('%i analyses failed out of %i\n',Nfails,Nanalyses)
failRate = Nfails/Nanalyses*100;
fprintf('This is a failure rate of %.2f%%\n',failRate)
NnegEigs = numel(find(eigsMatrix < 0.0));
fprintf('%i analyses with negative eigenvalues\n',NnegEigs)

uniqueEigs = unique(eigs);

% BY GM RECORD

failsEQ = zeros(1,NgMotions);
failsEQper = zeros(1,NgMotions);
negEigsEQ = zeros(1,NgMotions);
negEigsEQper = zeros(1,NgMotions);
runsPerEQ = Nanalyses/NgMotions;

for ii = 1:NgMotions
    
    failsEQ(ii) = sum(sum(sum(okMatrix(ii,:,:,:))))/-3;
    failsEQper(ii) = failsEQ(ii)/runsPerEQ*100;
    negEigsEQ(ii) = numel(find(eigsMatrix(ii,:,:,:) < 0.0));
    negEigsEQper(ii) = negEigsEQ(ii)/runsPerEQ*100;
    
end

% figure
% plot(1:NgMotions,failsEQper,'o--',1:NgMotions,negEigsEQper,'x--')
% grid on
% xlabel('EQ number')
% ylabel('Unstable analyses (%)')
% legend('Failed analysis','Negative eigenvalue','Location','Best')
% title('Instability by GM record')

% BY ASPECT RATIO

failsAR = zeros(1,NARs);
failsARper = zeros(1,NARs);
negEigsAR = zeros(1,NARs);
negEigsARper = zeros(1,NARs);
runsPerAR = Nanalyses/NARs;

for ii = 1:NARs
    
    failsAR(ii) = sum(sum(sum(okMatrix(:,ii,:,:))))/-3;
    failsARper(ii) = failsAR(ii)/runsPerAR*100;
    negEigsAR(ii) = numel(find(eigsMatrix(:,ii,:,:) < 0.0));
    negEigsARper(ii) = negEigsAR(ii)/runsPerEQ*100;
    
end

% figure
% plot(ARs,failsARper,'o--',ARs,negEigsARper,'x--')
% grid on
% xlabel('Aspect ratio')
% ylabel('Unstable analyses (%)')
% legend('Failed analysis','Negative eigenvalue','Location','Best')
% title('Instability by aspect ratio')

% BY PITCH

failsPitch = zeros(1,NPitchs);
failsPitchper = zeros(1,NPitchs);
negEigsPitch = zeros(1,NPitchs);
negEigsPitchper = zeros(1,NPitchs);
runsPerPitch = Nanalyses/NPitchs;

for ii = 1:NPitchs
    
    failsPitch(ii) = sum(sum(sum(okMatrix(:,:,:,ii))))/-3;
    failsPitchper(ii) = failsPitch(ii)/runsPerPitch*100;
    negEigsPitch(ii) = numel(find(eigsMatrix(:,:,:,ii) < 0.0));
    negEigsPitchper(ii) = negEigsPitch(ii)/runsPerEQ*100;
    
end

% figure
% plot(Pitchs*1000,failsPitchper,'o--',Pitchs*1000,negEigsPitchper,'x--')
% grid on
% xlabel('Pitch (mm)')
% ylabel('Unstable analyses (%)')
% legend('Failed analysis','Negative eigenvalue','Location','Best')
% title('Instability by pitch')

%% Some useful plots

% AVERAGE RESULTS OVER ALL GM RECORDS

% average demands over all GM records
avgdemandLHSMatrix = squeeze(mean(demandLHSMatrix,'omitnan'));
avgdemandRHSMatrix = squeeze(mean(demandRHSMatrix,'omitnan'));
avgratsLHSMatrix = squeeze(mean(ratsLHSMatrix,'omitnan'));
avgratsRHSMatrix = squeeze(mean(ratsRHSMatrix,'omitnan'));

avgdempersLHSMatrix = squeeze(mean(dempersLHSMatrix,'omitnan'));
avgdempersRHSMatrix = squeeze(mean(dempersRHSMatrix,'omitnan'));

% average displacements, acclerations and reactions over all GM records
avgdispxMatrix = squeeze(mean(dispxMatrix,'omitnan'));
avgdispyMatrix = squeeze(mean(dispyMatrix,'omitnan'));
avgdispzMatrix = squeeze(mean(dispzMatrix,'omitnan'));
avgaccelxMatrix = squeeze(mean(accelxMatrix,'omitnan'));
avgaccelyMatrix = squeeze(mean(accelyMatrix,'omitnan'));
avgreactxMatrix = squeeze(mean(reactxMatrix,'omitnan'));
avgreactyMatrix = squeeze(mean(reactyMatrix,'omitnan'));

% GEOMETRIC MEANS

% MAY NEED TO BE A LOOP THAT CONSIDERS NANS, AS THIS WILL AFFECT N VALUE
% AND PROD DOESN'T HAVE THE BUILT IN CAPABILITY OF OMITNAN

demandLHSMatrix2(demandLHSMatrix2==0) = 1;
demandRHSMatrix2(demandRHSMatrix2==0) = 1;

ratsLHSMatrix2(ratsLHSMatrix2==0) = 1;
ratsRHSMatrix2(ratsRHSMatrix2==0) = 1;

NMatrix(NMatrix==0) = 1;

geomeandemandLHSMatrix = nthroot(squeeze(prod(demandLHSMatrix2)),NMatrix);
geomeandemandRHSMatrix = nthroot(squeeze(prod(demandRHSMatrix2)),NMatrix);

geomeanratsLHSMatrix = nthroot(squeeze(prod(ratsLHSMatrix2)),NMatrix);
geomeanratsRHSMatrix = nthroot(squeeze(prod(ratsRHSMatrix2)),NMatrix);

geomeandempersLHSMatrix = nthroot(squeeze(prod(dempersLHSMatrix2)),NMatrix);
geomeandempersRHSMatrix = nthroot(squeeze(prod(dempersRHSMatrix2)),NMatrix);

proddispxMatrix2 = squeeze(prod(dispxMatrix2));
proddispyMatrix2 = squeeze(prod(dispyMatrix2));
proddispzMatrix2 = squeeze(prod(dispzMatrix2));
prodaccelxMatrix2 = squeeze(prod(accelxMatrix2));
prodaccelyMatrix2 = squeeze(prod(accelyMatrix2));
prodreactxMatrix2 = squeeze(prod(reactxMatrix2));
prodreactyMatrix2 = squeeze(prod(reactyMatrix2));

geomeandispxMatrix = zeros(Ndispx,NARs,NTns2,NPitchs);
geomeandispyMatrix = zeros(Ndispx,NARs,NTns2,NPitchs);
geomeandispzMatrix = zeros(Ndispx,NARs,NTns2,NPitchs);

for ii = 1:Ndispx
    
    geomeandispxMatrix(ii,:,:,:) = nthroot(squeeze(proddispxMatrix2(ii,:,:,:)),NMatrix);
    geomeandispyMatrix(ii,:,:,:) = nthroot(squeeze(proddispyMatrix2(ii,:,:,:)),NMatrix);
    geomeandispzMatrix(ii,:,:,:) = nthroot(squeeze(proddispzMatrix2(ii,:,:,:)),NMatrix);
    
end

geomeanaccelxMatrix = zeros(Naccelx,NARs,NTns2,NPitchs);
geomeanaccelyMatrix = zeros(Naccelx,NARs,NTns2,NPitchs);

for ii = 1:Naccelx
    
    geomeanaccelxMatrix(ii,:,:,:) = nthroot(squeeze(prodaccelxMatrix2(ii,:,:,:)),NMatrix);
    geomeanaccelyMatrix(ii,:,:,:) = nthroot(squeeze(prodaccelyMatrix2(ii,:,:,:)),NMatrix);
    
end

geomeanreactxMatrix = zeros(Nreactx,NARs,NTns2,NPitchs);
geomeanreactyMatrix = zeros(Nreactx,NARs,NTns2,NPitchs);

for ii = 1:Nreactx
    
    geomeanreactxMatrix(ii,:,:,:) = nthroot(squeeze(prodreactxMatrix2(ii,:,:,:)),NMatrix);
    geomeanreactyMatrix(ii,:,:,:) = nthroot(squeeze(prodreactyMatrix2(ii,:,:,:)),NMatrix);
    
end

defFromRot = squeeze(geomeandispzMatrix(1,:,:,:));%.*HMatrix;

% SELECT FIXED VARIABLES

% select fixed variable indices (too many variables to plot!)
AR2plot = 2;
pitch2plot = 3;
dispNode = 6; % -node 2 3 4 10 -dof 1 2 disp OR -node 1 2 3 4 9 10 -dof 1 2 3 disp
leftEdgeNode = 3;
rightEdgeNode = 4;
leftGNode = 5; % -node 1 2 3 4 5 6 7 8 9 10 -dof 1 2 reaction
rightGNode = 6;

driftMatrix = squeeze(dispxMatrix(:,dispNode,:,:,:))./HMatrixAll;
geomeandriftMatrix = squeeze(geomeandispxMatrix(dispNode,:,:,:))./HMatrix;

% True drift values (minus movement at base)
topDMatrix = squeeze(dispxMatrix(:,dispNode,:,:,:) - dispxMatrix(:,1,:,:,:));
topDMatrix(topDMatrix == 0) = 1;
geomeantopDMatrix = nthroot(squeeze(prod(topDMatrix)),NMatrix);

geoDrift = 100*squeeze(geomeandispxMatrix(dispNode,:,:,:))./HMatrix;
geoDrift2 = 100*geomeantopDMatrix./HMatrix;
geoRot = 100*defFromRot;

rackDemandLHSMatrix = ratsLHSMatrix.*PitchMatrix;
rackDemandRHSMatrix = ratsRHSMatrix.*PitchMatrix;

%% Geometric mean peak drift, total and from rotation (single AR, multiple pitches)

figure
h = plot(TnsCell{AR2plot},[geoDrift(AR2plot,1:NTnsArray(AR2plot),1); geoRot(AR2plot,1:NTnsArray(AR2plot),1)],'b--');
set(h, {'marker'}, {'o'; 's'});
hold on
h2 = plot(TnsCell{AR2plot},[geoDrift(AR2plot,1:NTnsArray(AR2plot),2); geoRot(AR2plot,1:NTnsArray(AR2plot),2)],'g--');
set(h2, {'marker'}, {'o'; 's'});
h3 = plot(TnsCell{AR2plot},[geoDrift(AR2plot,1:NTnsArray(AR2plot),3); geoRot(AR2plot,1:NTnsArray(AR2plot),3)],'r--');
set(h3, {'marker'}, {'o'; 's'});
h4 = plot(TnsCell{AR2plot},[geoDrift(AR2plot,1:NTnsArray(AR2plot),4); geoRot(AR2plot,1:NTnsArray(AR2plot),4)],'c--');
set(h4, {'marker'}, {'o'; 's'});
h5 = plot(TnsCell{AR2plot},[geoDrift(AR2plot,1:NTnsArray(AR2plot),5); geoRot(AR2plot,1:NTnsArray(AR2plot),5)],'k--');
set(h5, {'marker'}, {'o'; 's'});
grid on
ylim([0 inf])
xlabel('Period (s)')
ylabel('Drift (%)')
legend([h(1) h(2) h2(1) h3(1) h4(1) h5(1)],'Drift (total), Pitch = 1 mm',...
    'Drift (rotation)','Pitch = 2 mm','Pitch = 5 mm','Pitch = 10 mm',...
    'Pitch = 20 mm','Location','Best')
dim = [.2 .1 .3 .3];
ARno = num2str(ARs(AR2plot));
str = strcat({'Aspect ratio:'},{' '},{ARno});
annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',12,'FontWeight','bold');
fileLoc = strcat(filedir,'Figures/Drift_AR4');
saveas(gcf,fileLoc,'meta')
saveas(gcf,fileLoc,'fig')

%% Geometric mean peak drift, total (single pitch, multiple ARs)

figure
hold on
for ii = 1:NARs
plot(TnsCell{ii},100*squeeze(geomeantopDMatrix(ii,1:NTnsArray(ii),pitch2plot))./squeeze(HMatrix(ii,1:NTnsArray(ii),pitch2plot)),'o--')
end
grid on
ylim([0 inf])
xlabel('Period (s)')
ylabel('Drift (%)')
% title({'Geometric mean peak drift',['PT = ',num2str(PTper),'%']})
legend('AR = 2','AR = 4','AR = 6','AR = 8','Location','Best')
dim = [.2 .1 .3 .3];
Pitchno = num2str(1000*Pitchs(pitch2plot));
str = strcat({'Pitch:'},{' '},{Pitchno},{' mm'});
annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',12,'FontWeight','bold');

%% Geometric mean peak drift, total and from rotation (single pitch, multiple ARs)

myCols = ['b' 'r' 'g' 'm'];
figure
hold on
for ii = 1:NARs
h = plot(TnsCell{ii},[geoDrift2(ii,1:NTnsArray(ii),pitch2plot)' geoRot(ii,1:NTnsArray(ii),pitch2plot)'],'--');
set(h, {'marker'}, {'o'; 's'});
set(h, {'color'},{myCols(ii)});
end
grid on
ylim([0 inf])
xlabel('Period (s)')
ylabel('Drift (%)')
legend('AR = 2 (total)','AR = 2 (rotation)','AR = 4 (total)','AR = 4 (rotation)',...
    'AR = 6 (total)','AR = 6 (rotation)','AR = 8 (total)','AR = 8 (rotation)',...
    'Location','Best')
dim = [.7 .05 .3 .3];
Pitchno = num2str(1000*Pitchs(pitch2plot));
str = strcat({'Pitch:'},{' '},{Pitchno},{' mm'});
annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',12,'FontWeight','bold');
fileLoc = strcat(filedir,'Figures/Drift_P5');
saveas(gcf,fileLoc,'meta')
saveas(gcf,fileLoc,'fig')

proDriftFromRot = 100*squeeze(geoRot(:,:,pitch2plot))./squeeze(geoDrift2(:,:,pitch2plot));
disp(proDriftFromRot)

%% Rotation angle (single pitch, multiple ARs)

pitch2plot = 3;
figure
hold on
for ii = 1:NARs
plot(TnsCell{ii},defFromRot(ii,1:NTnsArray(ii),pitch2plot).*(180/pi),'o--')
end
grid on
ylim([0 inf])
xlabel('Period (s)')
ylabel('Rotation (degrees)')
% title({'Geometric mean peak drift',['PT = ',num2str(PTper),'%']})
legend('AR = 2','AR = 4','AR = 6','AR = 8','Location','Best')
dim = [.2 .1 .3 .3];
Pitchno = num2str(1000*Pitchs(pitch2plot));
str = strcat({'Pitch:'},{' '},{Pitchno},{' mm'});
annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',12,'FontWeight','bold');
fileLoc = strcat(filedir,'Figures/Rotation_P5');
saveas(gcf,fileLoc,'meta')
saveas(gcf,fileLoc,'fig')

%% Rotation angle (single AR, multiple pitches)

figure
plot(TnsCell{AR2plot},squeeze(defFromRot(AR2plot,1:NTnsArray(AR2plot),:)).*(180/pi),'o--')
grid on
ylim([0 inf])
xlabel('Period (s)')
ylabel('Rotation (degrees)')
legend('Pitch = 1 mm','Pitch = 2 mm','Pitch = 5 mm','Pitch = 10 mm',...
    'Pitch = 20 mm','Location','Best')
dim = [.2 .1 .3 .3];
ARno = num2str(ARs(AR2plot));
str = strcat({'Aspect ratio:'},{' '},{ARno});
annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',12,'FontWeight','bold');
fileLoc = strcat(filedir,'Figures/Rotation_AR4');
saveas(gcf,fileLoc,'meta')
saveas(gcf,fileLoc,'fig')

%% Geometric mean peak uplift, LHS (single AR, multiple pitches)
figure
plot(TnsCell{AR2plot},squeeze(geomeandispyMatrix(leftEdgeNode,AR2plot,1:NTnsArray(AR2plot),:)),'o--')
grid on
ylim([0 inf])
xlabel('Period (s)')
ylabel('Peak uplift LHS (m)')
legend('Pitch = 1 mm','Pitch = 2 mm','Pitch = 5 mm','Pitch = 10 mm',...
    'Pitch = 20 mm','Location','Best')
dim = [.2 .1 .3 .3];
ARno = num2str(ARs(AR2plot));
str = strcat({'Aspect ratio:'},{' '},{ARno});
annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',12,'FontWeight','bold');
fileLoc = strcat(filedir,'Figures/UpliftLHS_AR4');
saveas(gcf,fileLoc,'meta')
saveas(gcf,fileLoc,'fig')

%% Geometric mean peak uplift, RHS (single AR, multiple pitches)

figure
plot(TnsCell{AR2plot},squeeze(geomeandispyMatrix(rightEdgeNode,AR2plot,1:NTnsArray(AR2plot),:)),'o--')
grid on
ylim([0 inf])
xlabel('Period (s)')
ylabel('Peak uplift RHS (m)')
legend('Pitch = 1 mm','Pitch = 2 mm','Pitch = 5 mm','Pitch = 10 mm',...
    'Pitch = 20 mm','Location','Best')
dim = [.2 .1 .3 .3];
ARno = num2str(ARs(AR2plot));
str = strcat({'Aspect ratio:'},{' '},{ARno});
annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',12,'FontWeight','bold');
fileLoc = strcat(filedir,'Figures/UpliftRHS_AR4');
saveas(gcf,fileLoc,'meta')
saveas(gcf,fileLoc,'fig')

%% Geometric mean peak uplift, LHS (single pitch, multiple ARs)

figure
hold on
for ii = 1:NARs
    plot(TnsCell{ii},squeeze(geomeandispyMatrix(leftEdgeNode,ii,1:NTnsArray(ii),pitch2plot)),'o--')
end
grid on
ylim([0 inf])
xlabel('Period (s)')
ylabel('Peak uplift LHS (m)')
legend('AR = 2','AR = 4','AR = 6','AR = 8','Location','Best')
dim = [.2 .1 .3 .3];
Pitchno = num2str(1000*Pitchs(pitch2plot));
str = strcat({'Pitch:'},{' '},{Pitchno},{' mm'});
annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',12,'FontWeight','bold');
fileLoc = strcat(filedir,'Figures/UpliftLHS_P5');
saveas(gcf,fileLoc,'meta')
saveas(gcf,fileLoc,'fig')

%% Geometric mean peak uplift, RHS (single pitch, multiple ARs)

figure
hold on
for ii = 1:NARs
    plot(TnsCell{ii},squeeze(geomeandispyMatrix(rightEdgeNode,ii,1:NTnsArray(ii),pitch2plot)),'o--')
end
grid on
ylim([0 inf])
xlabel('Period (s)')
ylabel('Peak uplift RHS (m)')
legend('AR = 2','AR = 4','AR = 6','AR = 8','Location','Best')
dim = [.2 .1 .3 .3];
Pitchno = num2str(1000*Pitchs(pitch2plot));
str = strcat({'Pitch:'},{' '},{Pitchno},{' mm'});
annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',12,'FontWeight','bold');
fileLoc = strcat(filedir,'Figures/UpliftRHS_P5');
saveas(gcf,fileLoc,'meta')
saveas(gcf,fileLoc,'fig')

%% Geometric mean ratchet count, LHS (single AR, multiple pitches)

figure
plot(TnsCell{AR2plot},squeeze(geomeanratsLHSMatrix(AR2plot,1:NTnsArray(AR2plot),:)),'o--')
grid on
ylim([0 inf])
xlabel('Period (s)')
ylabel('No. of ratcheting actions LHS')
% title({'Geometric mean ratchet count (LHS)',['PT = ',num2str(PTper),'%']})
legend('Pitch = 1 mm','Pitch = 2 mm','Pitch = 5 mm','Pitch = 10 mm',...
    'Pitch = 20 mm','Location','Best')
dim = [.2 .1 .3 .3];
ARno = num2str(ARs(AR2plot));
str = strcat({'Aspect ratio:'},{' '},{ARno});
annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',12,'FontWeight','bold');
fileLoc = strcat(filedir,'Figures/RatchetCountLHS_AR4');
saveas(gcf,fileLoc,'meta')
saveas(gcf,fileLoc,'fig')

%% Geometric mean ratchet count, RHS (single AR, multiple pitches)

figure
plot(TnsCell{AR2plot},squeeze(geomeanratsRHSMatrix(AR2plot,1:NTnsArray(AR2plot),:)),'o--')
grid on
ylim([0 inf])
xlabel('Period (s)')
ylabel('No. of ratcheting actions RHS')
% title({'Geometric mean ratchet count (RHS)',['PT = ',num2str(PTper),'%']})
legend('Pitch = 1 mm','Pitch = 2 mm','Pitch = 5 mm','Pitch = 10 mm',...
    'Pitch = 20 mm','Location','Best')
dim = [.2 .1 .3 .3];
ARno = num2str(ARs(AR2plot));
str = strcat({'Aspect ratio:'},{' '},{ARno});
annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',12,'FontWeight','bold');
fileLoc = strcat(filedir,'Figures/RatchetCountRHS_AR4');
saveas(gcf,fileLoc,'meta')
saveas(gcf,fileLoc,'fig')

%% Geometric mean ratchet count, LHS (single pitch, multiple ARs)
figure
hold on
for ii = 1:NARs
    plot(TnsCell{ii},squeeze(geomeanratsLHSMatrix(ii,1:NTnsArray(ii),pitch2plot)),'o--')
end
grid on
ylim([0 inf])
xlabel('Period (s)')
ylabel('No. of ratcheting actions LHS')
legend('AR = 2','AR = 4','AR = 6','AR = 8','Location','Best')
dim = [.2 .1 .3 .3];
Pitchno = num2str(1000*Pitchs(pitch2plot));
str = strcat({'Pitch:'},{' '},{Pitchno},{' mm'});
annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',12,'FontWeight','bold');
fileLoc = strcat(filedir,'Figures/RatchetCountLHS_P5');
saveas(gcf,fileLoc,'meta')
saveas(gcf,fileLoc,'fig')

%% Geometric mean ratchet count, RHS (single pitch, multiple ARs)
figure
hold on
for ii = 1:NARs
    plot(TnsCell{ii},squeeze(geomeanratsRHSMatrix(ii,1:NTnsArray(ii),pitch2plot)),'o--')
end
grid on
ylim([0 inf])
xlabel('Period (s)')
ylabel('No. of ratcheting actions RHS')
legend('AR = 2','AR = 4','AR = 6','AR = 8','Location','Best')
dim = [.2 .1 .3 .3];
Pitchno = num2str(1000*Pitchs(pitch2plot));
str = strcat({'Pitch:'},{' '},{Pitchno},{' mm'});
annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',12,'FontWeight','bold');
fileLoc = strcat(filedir,'Figures/RatchetCountRHS_P5');
saveas(gcf,fileLoc,'meta')
saveas(gcf,fileLoc,'fig')

%% Geometric mean rack demand, LHS (single AR, multiple pitches)

figure
plot(TnsCell{AR2plot},squeeze(geomeanratsLHSMatrix(AR2plot,1:NTnsArray(AR2plot),:)).*(ones(NTnsArray(AR2plot),1)*Pitchs'),'o--')
grid on
ylim([0 inf])
xlabel('Period (s)')
ylabel('Rack demand LHS (m)')
% title({'Geometric mean rack demand (LHS)',['PT = ',num2str(PTper),'%']})
legend('Pitch = 1 mm','Pitch = 2 mm','Pitch = 5 mm','Pitch = 10 mm',...
    'Pitch = 20 mm','Location','Best')
dim = [.2 .1 .3 .3];
ARno = num2str(ARs(AR2plot));
str = strcat({'Aspect ratio:'},{' '},{ARno});
annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',12,'FontWeight','bold');
fileLoc = strcat(filedir,'Figures/RackDemandLHS_AR4');
saveas(gcf,fileLoc,'meta')
saveas(gcf,fileLoc,'fig')

%% Geometric mean rack demand, RHS (single AR, multiple pitches)

figure
plot(TnsCell{AR2plot},squeeze(geomeanratsRHSMatrix(AR2plot,1:NTnsArray(AR2plot),:)).*(ones(NTnsArray(AR2plot),1)*Pitchs'),'o--')
grid on
ylim([0 inf])
xlabel('Period (s)')
ylabel('Rack demand RHS (m)')
% title({'Geometric mean rack demand (RHS)',['PT = ',num2str(PTper),'%']})
legend('Pitch = 1 mm','Pitch = 2 mm','Pitch = 5 mm','Pitch = 10 mm',...
    'Pitch = 20 mm','Location','Best')
dim = [.2 .1 .3 .3];
ARno = num2str(ARs(AR2plot));
str = strcat({'Aspect ratio:'},{' '},{ARno});
annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',12,'FontWeight','bold');
fileLoc = strcat(filedir,'Figures/RackDemandRHS_AR4');
saveas(gcf,fileLoc,'meta')
saveas(gcf,fileLoc,'fig')

%% Geometric mean rack demand, LHS (single pitch, multiple ARs)

figure
hold on
for ii = 1:NARs
    plot(TnsCell{ii},squeeze(geomeanratsLHSMatrix(ii,1:NTnsArray(ii),pitch2plot))*Pitchs(pitch2plot),'o--')
end
grid on
ylim([0 inf])
xlabel('Period (s)')
ylabel('Rack demand LHS (m)')
legend('AR = 2','AR = 4','AR = 6','AR = 8','Location','Best')
dim = [.2 .1 .3 .3];
Pitchno = num2str(1000*Pitchs(pitch2plot));
str = strcat({'Pitch:'},{' '},{Pitchno},{' mm'});
annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',12,'FontWeight','bold');
fileLoc = strcat(filedir,'Figures/RackDemandLHS_P5');
saveas(gcf,fileLoc,'meta')
saveas(gcf,fileLoc,'fig')

%% Geometric mean rack demand, RHS (single pitch, multiple ARs)

figure
hold on
for ii = 1:NARs
    plot(TnsCell{ii},squeeze(geomeanratsRHSMatrix(ii,1:NTnsArray(ii),pitch2plot))*Pitchs(pitch2plot),'o--')
end
grid on
ylim([0 inf])
xlabel('Period (s)')
ylabel('Rack demand RHS (m)')
legend('AR = 2','AR = 4','AR = 6','AR = 8','Location','Best')
dim = [.2 .1 .3 .3];
Pitchno = num2str(1000*Pitchs(pitch2plot));
str = strcat({'Pitch:'},{' '},{Pitchno},{' mm'});
annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',12,'FontWeight','bold');
fileLoc = strcat(filedir,'Figures/RackDemandRHS_P5');
saveas(gcf,fileLoc,'meta')
saveas(gcf,fileLoc,'fig')

%% Geometreic mean inelastic demand, LHS (single AR, multiple pitches)

figure
plot(TnsCell{AR2plot},squeeze(geomeandemandLHSMatrix(AR2plot,1:NTnsArray(AR2plot),:)),'o--')
grid on
ylim([0 inf])
xlabel('Period (s)')
ylabel('Inelastic displacement LHS (m)')
% title({'Geometric mean cumulative inelastc dissipater displacement (LHS)',['PT = ',num2str(PTper),'%']})
legend('Pitch = 1 mm','Pitch = 2 mm','Pitch = 5 mm','Pitch = 10 mm',...
    'Pitch = 20 mm','Location','Best')
dim = [.2 .1 .3 .3];
ARno = num2str(ARs(AR2plot));
str = strcat({'Aspect ratio:'},{' '},{ARno});
annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',12,'FontWeight','bold');
fileLoc = strcat(filedir,'Figures/DemandLHS_AR4');
saveas(gcf,fileLoc,'meta')
saveas(gcf,fileLoc,'fig')

%% Geometreic mean inelastic demand, RHS (single AR, multiple pitches)

figure
plot(TnsCell{AR2plot},squeeze(geomeandemandRHSMatrix(AR2plot,1:NTnsArray(AR2plot),:)),'o--')
grid on
ylim([0 inf])
xlabel('Period (s)')
ylabel('Inelastic displacement RHS (m)')
% title({'Geometric mean cumulative inelastc dissipater displacement (RHS)',['PT = ',num2str(PTper),'%']})
legend('Pitch = 1 mm','Pitch = 2 mm','Pitch = 5 mm','Pitch = 10 mm',...
    'Pitch = 20 mm','Location','Best')
dim = [.2 .1 .3 .3];
ARno = num2str(ARs(AR2plot));
str = strcat({'Aspect ratio:'},{' '},{ARno});
annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',12,'FontWeight','bold');
fileLoc = strcat(filedir,'Figures/DemandRHS_AR4');
saveas(gcf,fileLoc,'meta')
saveas(gcf,fileLoc,'fig')

%% Geometric mean inelastic demand LHS (single pitch, multiple ARs)

figure
hold on
for ii = 1:NARs
    plot(TnsCell{ii},squeeze(geomeandemandLHSMatrix(ii,1:NTnsArray(ii),pitch2plot)),'o--')
end
grid on
ylim([0 inf])
xlabel('Period (s)')
ylabel('Inelastic displacement LHS (m)')
legend('AR = 2','AR = 4','AR = 6','AR = 8','Location','Best')
dim = [.2 .1 .3 .3];
Pitchno = num2str(1000*Pitchs(pitch2plot));
str = strcat({'Pitch:'},{' '},{Pitchno},{' mm'});
annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',12,'FontWeight','bold');
fileLoc = strcat(filedir,'Figures/DemandLHS_P5');
saveas(gcf,fileLoc,'meta')
saveas(gcf,fileLoc,'fig')

%% Geometric mean inelastic demand RHS (single pitch, multiple ARs)

figure
hold on
for ii = 1:NARs
    plot(TnsCell{ii},squeeze(geomeandemandRHSMatrix(ii,1:NTnsArray(ii),pitch2plot)),'o--')
end
grid on
ylim([0 inf])
xlabel('Period (s)')
ylabel('Inelastic displacement RHS (m)')
legend('AR = 2','AR = 4','AR = 6','AR = 8','Location','Best')
dim = [.2 .1 .3 .3];
Pitchno = num2str(1000*Pitchs(pitch2plot));
str = strcat({'Pitch:'},{' '},{Pitchno},{' mm'});
annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',12,'FontWeight','bold');
fileLoc = strcat(filedir,'Figures/DemandRHS_P5');
saveas(gcf,fileLoc,'meta')
saveas(gcf,fileLoc,'fig')

%% Comparison of inelastic demand, LHS (single AR, multiple pitches)
% AR2plot = 4;
figure
plot(TnsCell{AR2plot},squeeze(geomeandempersLHSMatrix(AR2plot,1:NTnsArray(AR2plot),:)),'o--')
grid on
ylim([0 inf])
xlabel('Period (s)')
ylabel('Demand / Demand_P_ _=_ _1_m_m LHS')
% title({'Effect of pitch on demand (LHS)',['PT = ',num2str(PTper),'%']})
legend('Pitch = 1 mm','Pitch = 2 mm','Pitch = 5 mm','Pitch = 10 mm',...
    'Pitch = 20 mm','Location','Best')
dim = [.2 .1 .3 .3];
ARno = num2str(ARs(AR2plot));
str = strcat({'Aspect ratio:'},{' '},{ARno});
annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',12,'FontWeight','bold');
fileLoc = strcat(filedir,'Figures/RelativeDemandLHS_AR4');
saveas(gcf,fileLoc,'meta')
saveas(gcf,fileLoc,'fig')

%% Comparison of inelastic demand, RHS (single AR, multiple pitches)

figure
plot(TnsCell{AR2plot},squeeze(geomeandempersRHSMatrix(AR2plot,1:NTnsArray(AR2plot),:)),'o--')
grid on
ylim([0 inf])
xlabel('Period (s)')
ylabel('Demand / Demand_P_ _=_ _1_m_m RHS')
% title({'Effect of pitch on demand (RHS)',['PT = ',num2str(PTper),'%']})
legend('Pitch = 1 mm','Pitch = 2 mm','Pitch = 5 mm','Pitch = 10 mm',...
    'Pitch = 20 mm','Location','Best')
dim = [.2 .1 .3 .3];
ARno = num2str(ARs(AR2plot));
str = strcat({'Aspect ratio:'},{' '},{ARno});
annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',12,'FontWeight','bold');
fileLoc = strcat(filedir,'Figures/RelativeDemandRHS_AR4');
saveas(gcf,fileLoc,'meta')
saveas(gcf,fileLoc,'fig')

%% Comparison of inelastic demand (All)

for ii = 1:NARs
    
    figure
    plot(TnsCell{ii},squeeze(geomeandempersLHSMatrix(ii,1:NTnsArray(ii),:)),'o--')
    grid on
    ylim([0 inf])
    xlabel('Period (s)')
    ylabel('Demand / Demand_P_ _=_ _1_m_m LHS')
    legend('Pitch = 1 mm','Pitch = 2 mm','Pitch = 5 mm','Pitch = 10 mm',...
        'Pitch = 20 mm','Location','Best')
    dim = [.2 .1 .3 .3];
    ARno = num2str(ARs(ii));
    str = strcat({'Aspect ratio:'},{' '},{ARno});
    annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',12,'FontWeight','bold');
    
    figure
    plot(TnsCell{ii},squeeze(geomeandempersRHSMatrix(ii,1:NTnsArray(ii),:)),'o--')
    grid on
    ylim([0 inf])
    xlabel('Period (s)')
    ylabel('Demand / Demand_P_ _=_ _1_m_m RHS')
    legend('Pitch = 1 mm','Pitch = 2 mm','Pitch = 5 mm','Pitch = 10 mm',...
        'Pitch = 20 mm','Location','Best')
    dim = [.2 .1 .3 .3];
    ARno = num2str(ARs(ii));
    str = strcat({'Aspect ratio:'},{' '},{ARno});
    annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',12,'FontWeight','bold');
    
end

%% PEAK UPLIFT AND INELASTIC DEMAND

upliftLHS = squeeze(dispyMatrix(:,leftEdgeNode,:,:,:));
upliftRHS = squeeze(dispyMatrix(:,rightEdgeNode,:,:,:));

moneyShotLHS = demandLHSMatrix./upliftLHS;
moneyShotRHS = demandRHSMatrix./upliftRHS;

suspectCasesLHS = moneyShotLHS < 0.8;
suspectCasesRHS = moneyShotRHS < 0.8;

maxDriftAll = max(max(max(max(max(abs(squeeze(dispxMatrix(:,dispNode,:,:,:)))./HMatrixAll)))));
maxDriftAllPlan = max(max(max(max(max(abs(squeeze(dispxMatrix(:,dispNode,:,:,:)))/(2*B))))));
maxDemandLHS = max(abs(demandLHSMatrix(:)));
maxDemandRHS = max(abs(demandRHSMatrix(:)));

maxDemandLocation = find(demandLHSMatrix == maxDemandLHS);
maxDriftLocation = find((squeeze(dispxMatrix(:,dispNode,:,:,:))./HMatrixAll) == maxDriftAll);

% Convert location to indices!!
arraySizes = [NgMotions, NARs, NTns, NPitchs];
[I, J, K] = ind2sub(arraySizes,maxDriftLocation);

theDemand = demandLHSMatrix(I,J,K);

%% Inelastic dissipater demand vs. peak uplift, LHS

figure
plot(upliftLHS(:),demandLHSMatrix(:),'b.')
grid on
hold on
plot([0 0.5],[0 0.5],'g--')
plot([0 0.3],[0 3],'r--')
legend('Data','Ratio = 1','Ratio = 10','Location','NorthWest')
xlabel('Peak uplift LHS (m)')
ylabel('Inelastic dissipater demand LHS (m)')
% title('LHS dissipater demand ratios')
fileLoc = strcat(filedir,'Figures/DemandUpliftLHS');
saveas(gcf,fileLoc,'meta')
saveas(gcf,fileLoc,'fig')

%% Inelastic dissipater demand vs. peak uplift, RHS

figure
plot(upliftRHS(:),demandRHSMatrix(:),'b.')
grid on
hold on
plot([0 0.5],[0 0.5],'g--')
plot([0 0.3],[0 3],'r--')
legend('Data','Ratio = 1','Ratio = 10','Location','NorthWest')
xlabel('Peak uplift RHS (m)')
ylabel('Inelastic dissipater demand RHS (m)')
% title('RHS dissipater demand ratios')
fileLoc = strcat(filedir,'Figures/DemandUpliftRHS');
saveas(gcf,fileLoc,'meta')
saveas(gcf,fileLoc,'fig')

%% Inelastic dissipater demand vs. peak uplift stats
% Initial stats

ratioArrayLHS = [];
ratioArrayRHS = [];

for ii = 1:NARs
    
    forThisLHS = moneyShotLHS(:,ii,1:NTnsArray(ii),:);
    ratioArrayLHS = [ratioArrayLHS; forThisLHS(:)];
    
    forThisRHS = moneyShotRHS(:,ii,1:NTnsArray(ii),:);
    ratioArrayRHS = [ratioArrayRHS; forThisRHS(:)];
    
end

ratioArray = [ratioArrayLHS ratioArrayRHS];

demandStatsFunc = @(x) [min(x) max(x) mean(x) geomean(x) std(x)];

demandStatsLHS = demandStatsFunc(ratioArrayLHS);
demandStatsRHS = demandStatsFunc(ratioArrayRHS);
demandStats = [demandStatsLHS; demandStatsRHS];

% fprintf('Dissipater demand ratio statistics (LHS, RHS):\n')
% fprintf('Min: %.2f, %.2f\n',demandStats(:,1))
% fprintf('Max: %.2f, %.2f\n',demandStats(:,2))
% fprintf('Mean: %.2f, %.2f\n',demandStats(:,3))
% fprintf('Geometric mean: %.2f, %.2f\n',demandStats(:,4))
% fprintf('Standard deviation: %.2f, %.2f\n',demandStats(:,5))

%% Inelastic dissipater demand vs. peak uplift stats
% Normal distribution (Gaussian, or additive normal, distribution)

% Omit zero results (not required for normal distribution analysis)
omitZerosLHS = ratioArrayLHS;%(ratioArrayLHS ~= 0);
omitZerosRHS = ratioArrayRHS;%(ratioArrayRHS ~= 0);

% Sample size
NLHSnorm = length(omitZerosLHS);
NRHSnorm = length(omitZerosRHS);

% Arithmetic mean
meanLHS = mean(omitZerosLHS);
meanRHS = mean(omitZerosRHS);

% Additive standard deviation
stdLHS = std(omitZerosLHS);
stdRHS = std(omitZerosRHS);

% Interval of confidence
% 68.3%
IC1LHSnorm = [meanLHS-stdLHS meanLHS+stdLHS];
IC1RHSnorm = [meanRHS-stdRHS meanRHS+stdRHS];
% 95.5%
IC2LHSnorm = [meanLHS-(2*stdLHS) meanLHS+(2*stdLHS)];
IC2RHSnorm = [meanRHS-(2*stdRHS) meanRHS+(2*stdRHS)];
% 99.7%
IC3LHSnorm = [meanLHS-(3*stdLHS) meanLHS+(3*stdLHS)];
IC3RHSnorm = [meanRHS-(3*stdRHS) meanRHS+(3*stdRHS)];

% True IC representation
repIC1LHSnorm = sum((omitZerosLHS < IC1LHSnorm(2)) & (omitZerosLHS > IC1LHSnorm(1)))/NLHSnorm*100;
repIC1RHSnorm = sum((omitZerosRHS < IC1RHSnorm(2)) & (omitZerosRHS > IC1RHSnorm(1)))/NRHSnorm*100;
repIC2LHSnorm = sum((omitZerosLHS < IC2LHSnorm(2)) & (omitZerosLHS > IC2LHSnorm(1)))/NLHSnorm*100;
repIC2RHSnorm = sum((omitZerosRHS < IC2RHSnorm(2)) & (omitZerosRHS > IC2RHSnorm(1)))/NRHSnorm*100;
repIC3LHSnorm = sum((omitZerosLHS < IC3LHSnorm(2)) & (omitZerosLHS > IC3LHSnorm(1)))/NLHSnorm*100;
repIC3RHSnorm = sum((omitZerosRHS < IC3RHSnorm(2)) & (omitZerosRHS > IC3RHSnorm(1)))/NRHSnorm*100;
repICnorm = [repIC1LHSnorm repIC1RHSnorm;repIC2RHSnorm repIC2RHSnorm;repIC3RHSnorm repIC3RHSnorm];

% repICnorm = [repIC1LHSnorm repIC1RHSnorm;repIC2LHSnorm repIC2RHSnorm;repIC3LHSnorm repIC3RHSnorm];
% fprintf('Confidence intervals using normal distribution properties:\n')
% disp(repICnorm)

% Portion of sample less than upper bound (captured by size recommendation)
lessThan1Stdnorm = [sum(omitZerosLHS < IC1LHSnorm(2))/NLHSnorm*100 sum(omitZerosRHS < IC1RHSnorm(2))/NRHSnorm*100];
lessThan2Stdnorm = [sum(omitZerosLHS < IC2LHSnorm(2))/NLHSnorm*100 sum(omitZerosRHS < IC2RHSnorm(2))/NRHSnorm*100];
lessThan3Stdnorm = [sum(omitZerosLHS < IC3LHSnorm(2))/NLHSnorm*100 sum(omitZerosRHS < IC3RHSnorm(2))/NRHSnorm*100];

% lessThanStdnorm = [lessThan1Stdnorm; lessThan2Stdnorm; lessThan3Stdnorm];
% fprintf('Demand cases captured by normal distribution properties:\n')
% disp(lessThanStdnorm)

%% Inelastic dissipater demand vs. peak uplift stats
% Log-normal distribution (multiplicative normal distribution)

% Omit zero results (required for log-normal distribution analysis)
omitZerosLHS = ratioArrayLHS(ratioArrayLHS ~= 0);
omitZerosRHS = ratioArrayRHS(ratioArrayRHS ~= 0);

% Sample size
NLHS = length(omitZerosLHS);
NRHS = length(omitZerosRHS);

% Geometric mean
geomeanRatioLHS_ = geomean(omitZerosLHS);
geomeanRatioRHS_ = geomean(omitZerosRHS);

% Multiplicative standard deviation
multiStdLHS = exp(sqrt((sum((log(omitZerosLHS/geomeanRatioLHS_)).^2))/(NLHS-1)));
multiStdRHS = exp(sqrt((sum((log(omitZerosRHS/geomeanRatioRHS_)).^2))/(NRHS-1)));

% Interval of confidence
% 68.3%
IC1LHS = [geomeanRatioLHS_/multiStdLHS geomeanRatioLHS_*multiStdLHS];
IC1RHS = [geomeanRatioRHS_/multiStdRHS geomeanRatioRHS_*multiStdRHS];
% 95.5%
IC2LHS = [geomeanRatioLHS_/(multiStdLHS^2) geomeanRatioLHS_*(multiStdLHS^2)];
IC2RHS = [geomeanRatioRHS_/(multiStdRHS^2) geomeanRatioRHS_*(multiStdRHS^2)];
% 99.7%
IC3LHS = [geomeanRatioLHS_/(multiStdLHS^3) geomeanRatioLHS_*(multiStdLHS^3)];
IC3RHS = [geomeanRatioRHS_/(multiStdRHS^3) geomeanRatioRHS_*(multiStdRHS^3)];

% Display summary statistics
fprintf('Dissipater demand ratio statistics (LHS, RHS):\n')
fprintf('Min: %.2f, %.2f\n',demandStats(:,1))
fprintf('Max: %.2f, %.2f\n',demandStats(:,2))
fprintf('Geometric mean: %.2f, %.2f\n',geomeanRatioLHS_,geomeanRatioRHS_)
fprintf('Multiplicative standard deviation: %.2f, %.2f\n',multiStdLHS,multiStdRHS)

% Display distribution limits
fprintf('Log-normal distribution limits of rack demand:\n')
fprintf('One standard deviation:\n')
disp([IC1LHS; IC1RHS])
fprintf('Two standard deviations:\n')
disp([IC2LHS; IC2RHS])
fprintf('Three standard deviations:\n')
disp([IC3LHS; IC3RHS])

% True IC representation
repIC1LHS = sum((omitZerosLHS < IC1LHS(2)) & (omitZerosLHS > IC1LHS(1)))/NLHS*100;
repIC1RHS = sum((omitZerosRHS < IC1RHS(2)) & (omitZerosRHS > IC1RHS(1)))/NRHS*100;
repIC2LHS = sum((omitZerosLHS < IC2LHS(2)) & (omitZerosLHS > IC2LHS(1)))/NLHS*100;
repIC2RHS = sum((omitZerosRHS < IC2RHS(2)) & (omitZerosRHS > IC2RHS(1)))/NRHS*100;
repIC3LHS = sum((omitZerosLHS < IC3LHS(2)) & (omitZerosLHS > IC3LHS(1)))/NLHS*100;
repIC3RHS = sum((omitZerosRHS < IC3RHS(2)) & (omitZerosRHS > IC3RHS(1)))/NRHS*100;

repIC = [repIC1LHS repIC1RHS;repIC2LHS repIC2RHS;repIC3LHS repIC3RHS];
fprintf('Confidence intervals using log-normal distribution properties:\n')
disp(repIC)

% Portion of sample less than upper bound (captured by size recommendation)
lessThan1Std = [sum(omitZerosLHS < IC1LHS(2))/NLHS*100 sum(omitZerosRHS < IC1RHS(2))/NRHS*100];
lessThan2Std = [sum(omitZerosLHS < IC2LHS(2))/NLHS*100 sum(omitZerosRHS < IC2RHS(2))/NRHS*100];
lessThan3Std = [sum(omitZerosLHS < IC3LHS(2))/NLHS*100 sum(omitZerosRHS < IC3RHS(2))/NRHS*100];

lessThanStd = [lessThan1Std; lessThan2Std; lessThan3Std];
fprintf('Demand cases captured by log-normal distribution properties:\n')
disp(lessThanStd)

%% Inelastic dissipater demand vs. peak uplift stats
% Log-normal distribution (multiplicative normal distribution) by pitch

geomeanRatioArrayLHS = zeros(1,NPitchs);
geomeanRatioArrayRHS = zeros(1,NPitchs);

multiStdArrayLHS = zeros(1,NPitchs);
multiStdArrayRHS = zeros(1,NPitchs);

IC1LHS_ = zeros(2,NPitchs); 
IC1RHS_ = zeros(2,NPitchs);
IC2LHS_ = zeros(2,NPitchs);
IC2RHS_ = zeros(2,NPitchs);
IC3LHS_ = zeros(2,NPitchs);
IC3RHS_ = zeros(2,NPitchs);

for ii = 1:NPitchs
    
    andNowThisLHS = [];
    andNowThisRHS = [];
    
    for jj = 1:NARs
        
        andNowThisLHSInner = squeeze(moneyShotLHS(:,jj,1:NTnsArray(jj),ii));
        andNowThisLHS = [andNowThisLHS andNowThisLHSInner(:)'];
        
        andNowThisRHSInner = squeeze(moneyShotRHS(:,jj,1:NTnsArray(jj),ii));
        andNowThisRHS = [andNowThisRHS andNowThisRHSInner(:)']; 
        
    end
    
    % Omit zero results (required for log-normal distribution analysis)
    omitZerosLHS = andNowThisLHS(andNowThisLHS ~= 0);
    omitZerosRHS = andNowThisRHS(andNowThisRHS ~= 0);
    
    % Sample size
    NLHS = length(omitZerosLHS);
    NRHS = length(omitZerosRHS);
    
    % Geometric mean
    geomeanRatioArrayLHS(ii) = geomean(omitZerosLHS);
    geomeanRatioArrayRHS(ii) = geomean(omitZerosRHS);
    
    % Multiplicative standard deviation
    multiStdArrayLHS(ii) = exp(sqrt((sum((log(omitZerosLHS/geomeanRatioArrayLHS(ii))).^2))/(NLHS-1)));
    multiStdArrayRHS(ii) = exp(sqrt((sum((log(omitZerosRHS/geomeanRatioArrayRHS(ii))).^2))/(NRHS-1)));
    
    % Interval of confidence
    % 68.3%
    IC1LHS_(:,ii) = [geomeanRatioArrayLHS(ii)/multiStdArrayLHS(ii) geomeanRatioArrayLHS(ii)*multiStdArrayLHS(ii)];
    IC1RHS_(:,ii) = [geomeanRatioArrayRHS(ii)/multiStdArrayRHS(ii) geomeanRatioArrayLHS(ii)*multiStdArrayRHS(ii)];
    % 95.5%
    IC2LHS_(:,ii) = [geomeanRatioArrayLHS(ii)/(multiStdArrayLHS(ii)^2) geomeanRatioArrayLHS(ii)*(multiStdArrayLHS(ii)^2)];
    IC2RHS_(:,ii) = [geomeanRatioArrayRHS(ii)/(multiStdArrayRHS(ii)^2) geomeanRatioArrayLHS(ii)*(multiStdArrayRHS(ii)^2)];
    % 99.7%
    IC3LHS_(:,ii) = [geomeanRatioArrayLHS(ii)/(multiStdArrayLHS(ii)^3) geomeanRatioArrayLHS(ii)*(multiStdArrayLHS(ii)^3)];
    IC3RHS_(:,ii) = [geomeanRatioArrayRHS(ii)/(multiStdArrayRHS(ii)^3) geomeanRatioArrayLHS(ii)*(multiStdArrayRHS(ii)^3)];
    
end

geomeanRatioArray = [geomeanRatioArrayLHS; geomeanRatioArrayRHS];
multiStdArray = [multiStdArrayLHS; multiStdArrayRHS];
relativeGeomeanRatioArrayLHS = 100*geomeanRatioArrayLHS./geomeanRatioArrayLHS(1);
relativeGeomeanRatioArrayRHS = 100*geomeanRatioArrayRHS./geomeanRatioArrayRHS(1);
relativeGeomeanRatioArray = [relativeGeomeanRatioArrayLHS; relativeGeomeanRatioArrayRHS];

% Dissipater demand summary statistics by pitch
fprintf('Pitch size influence on dissipater demand ratio:\n')
fprintf('Geometric mean:\n')
disp(geomeanRatioArray)
fprintf('Multiplicative Standard deviation:\n')
disp(multiStdArray)
fprintf('Geometric mean as a percentage of 1mm pitch value:\n')
disp(relativeGeomeanRatioArray)

%% Inelastic dissipater demand vs. peak uplift stats
% Demand ratios by pitch

NStats = size(demandStats,2);
StatsArrayLHS = zeros(NARs,NTns2,NPitchs,NStats);
StatsArrayRHS = zeros(NARs,NTns2,NPitchs,NStats);
StatsArray = zeros(2,NARs,NTns2,NPitchs,NStats);

ctr = 0;
for ii = 1:NARs
    
    for jj = 1:NTnsArray(ii)
        
        for kk = 1:NPitchs
            
            thisDataLHS = moneyShotLHS(:,ii,jj,kk);
            thisDataRHS = moneyShotRHS(:,ii,jj,kk);
            
            StatsArrayLHS(ii,jj,kk,:) = demandStatsFunc(thisDataLHS);
            StatsArrayRHS(ii,jj,kk,:) = demandStatsFunc(thisDataRHS);
            
            % Omit zero results
            omitZerosLHS = thisDataLHS(thisDataLHS ~= 0);
            omitZerosRHS = thisDataRHS(thisDataRHS ~= 0);
            
            % Geometric mean            
            geomeanLHS = geomean(omitZerosLHS);
            geomeanRHS = geomean(omitZerosRHS);
            
            StatsArrayLHS(ii,jj,kk,4) = geomeanLHS;
            StatsArrayRHS(ii,jj,kk,4) = geomeanRHS;
            
            StatsArray(:,ii,jj,kk,:) = [StatsArrayLHS(ii,jj,kk,:) StatsArrayRHS(ii,jj,kk,:)];
            
            ctr = ctr + 1;
        end
        
    end
    
end

%% Demand ratios by pitch (one structure per aspect ratio)

ARTn2plot = [2 4 5 7];
% ARTn2plot = [1 3 3 4];

figure
h = plot(Pitchs*1000,[squeeze(StatsArrayLHS(1,ARTn2plot(1),:,4)) squeeze(StatsArrayRHS(1,ARTn2plot(1),:,4))],'bo--');
set(h, {'marker'}, {'o'; 's'});
hold on
h2 = plot(Pitchs*1000,[squeeze(StatsArrayLHS(2,ARTn2plot(2),:,4)) squeeze(StatsArrayRHS(2,ARTn2plot(2),:,4))],'go--');
set(h2, {'marker'}, {'o'; 's'});
h3 = plot(Pitchs*1000,[squeeze(StatsArrayLHS(3,ARTn2plot(3),:,4)) squeeze(StatsArrayRHS(3,ARTn2plot(3),:,4))],'ro--');
set(h3, {'marker'}, {'o'; 's'});
h4 = plot(Pitchs*1000,[squeeze(StatsArrayLHS(4,ARTn2plot(4),:,4)) squeeze(StatsArrayRHS(4,ARTn2plot(4),:,4))],'ko--');
set(h4, {'marker'}, {'o'; 's'});
grid on
ylim([0 inf])
xlabel('Pitch (mm)')
ylabel('Geometric mean of demand ratio')

legData = [];
for ii = 1:NARs
    Tnno = num2str(TnsCell{ii}(ARTn2plot(ii)));
    myEntryLHS = strcat({'AR = '},{num2str(ARs(ii))},{', Tn = '},{Tnno},' s (LHS)');
    myEntryRHS = strcat({'AR = '},{num2str(ARs(ii))},{', Tn = '},{Tnno},' s (RHS)');
    legData = [legData myEntryLHS myEntryRHS];
end
legend(legData,'Location','Best')
fileLoc = strcat(filedir,'Figures/DemandRatio_T04_08_12_16');
saveas(gcf,fileLoc,'meta')
saveas(gcf,fileLoc,'fig')

% Showing all records
figure
h = plot(Pitchs*1000,squeeze(StatsArrayLHS(1,2,:,4)),'go--','LineWidth',3);
hold on
h2 = plot(Pitchs*1000,squeeze(moneyShotLHS(:,1,2,:)),'ko--','Color',[0 0 0]+0.6);
h = plot(Pitchs*1000,squeeze(StatsArrayLHS(1,2,:,4)),'go--','LineWidth',3);
grid on
ylim([0 inf])
xlabel('Pitch (mm)')
ylabel('Demand ratio')
legend('Geometric mean','Individual analyses','Location','Best')
dim = [.15 .08 .1 .1];
str = '   AR = 2, Tn = 0.4';
annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',12,'FontWeight','bold');

%% Demand ratio, LHS (All)

for ii = 1:NARs
    
    legData = {};
    figure
    hold on
    
    for jj = 1:NTnsArray(ii)
        h = plot(Pitchs*1000,squeeze(StatsArrayLHS(ii,jj,:,4)),'o--');
        Tnno = num2str(TnsCell{ii}(jj));
        myEntryLHS = strcat({'Tn = '},{Tnno},' s (LHS)');
        legData = [legData myEntryLHS];
    end
    
    dim = [.2 .2 .3 .3];
    ARno = num2str(ARs(ii));
    str = strcat({'Aspect ratio:'},{' '},{ARno});
    annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',12,'FontWeight','bold');
    
    grid on
    ylim([0 inf])
    xlabel('Pitch (mm)')
    ylabel('Geometric mean of demand ratio')
    legend(legData,'Location','Best')
    
    fileLoc = strcat(filedir,'Figures/DemandRatio_AR',num2str(ARs(ii)));
    saveas(gcf,fileLoc,'meta')
    saveas(gcf,fileLoc,'fig')
    
end

%% Rack demand ratio statistics
% Rack demand vs. peak uplift stats

rackDemandRatioLHSMatrix = rackDemandLHSMatrix./upliftLHS;
rackDemandRatioRHSMatrix = rackDemandRHSMatrix./upliftRHS;

rackDemandRatioArrayLHS = [];
rackDemandRatioArrayRHS = [];

for ii = 1:NARs
    
    forThisLHS = rackDemandRatioLHSMatrix(:,ii,1:NTnsArray(ii),:);
    rackDemandRatioArrayLHS = [rackDemandRatioArrayLHS; forThisLHS(:)];
    
    forThisRHS = rackDemandRatioRHSMatrix(:,ii,1:NTnsArray(ii),:);
    rackDemandRatioArrayRHS = [rackDemandRatioArrayRHS; forThisRHS(:)];
    
end

rackDemandRatioArray = [rackDemandRatioArrayLHS; rackDemandRatioArrayRHS];

rackDemandStats2 = demandStatsFunc(rackDemandRatioArray);

rackDemandStatsLHS = demandStatsFunc(rackDemandRatioArrayLHS);
rackDemandStatsRHS = demandStatsFunc(rackDemandRatioArrayRHS);
rackDemandStats = [rackDemandStatsLHS; rackDemandStatsRHS];

statsByPitchLHS = zeros(NPitchs,5);

for ii = 1:NPitchs
    
    andNowThis = [];
    
    for jj = 1:NARs
        
        andNowThisInner = squeeze(rackDemandRatioLHSMatrix(:,jj,1:NTnsArray(jj),ii));
        andNowThis = [andNowThis andNowThisInner(:)']; 
        
    end
    
    statsByPitchLHS(ii,:) = demandStatsFunc(andNowThis(:));
    
end

% fprintf('Rack demand ratio statistics (LHS, RHS):\n')
% fprintf('Min: %.2f, %.2f\n',rackDemandStats(:,1))
% fprintf('Max: %.2f, %.2f\n',rackDemandStats(:,2))
% fprintf('Mean: %.2f, %.2f\n',rackDemandStats(:,3))
% fprintf('Geometric mean: %.2f, %.2f\n',rackDemandStats(:,4))
% fprintf('Standard deviation: %.2f, %.2f\n',rackDemandStats(:,5))

%% Rack demand ratio statistics
% Log-normal distribution (multiplicative normal distribution)

% Omit zero results (required for log-normal distribution analysis)
omitZerosLHS = rackDemandRatioArrayLHS(rackDemandRatioArrayLHS ~= 0);
omitZerosRHS = rackDemandRatioArrayRHS(rackDemandRatioArrayRHS ~= 0);
omitZeros = [omitZerosLHS; omitZerosRHS];

% Sample size
NLHS = length(omitZerosLHS);
NRHS = length(omitZerosRHS);
NAll = length(omitZeros);

% Geometric mean
geomeanRackRatioLHS_ = geomean(omitZerosLHS);
geomeanRackRatioRHS_ = geomean(omitZerosRHS);
geomeanRackRatio_ = geomean(omitZeros);

% Multiplicative standard deviation
multiStdRackLHS = exp(sqrt((sum((log(omitZerosLHS/geomeanRackRatioLHS_)).^2))/(NLHS-1)));
multiStdRackRHS = exp(sqrt((sum((log(omitZerosRHS/geomeanRackRatioRHS_)).^2))/(NRHS-1)));
multiStdRack = exp(sqrt((sum((log(omitZeros/geomeanRackRatio_)).^2))/(NAll-1)));

% Interval of confidence
% 68.3%
IC1RackLHS = [geomeanRackRatioLHS_/multiStdRackLHS geomeanRackRatioLHS_*multiStdRackLHS];
IC1RackRHS = [geomeanRackRatioRHS_/multiStdRackRHS geomeanRackRatioRHS_*multiStdRackRHS];
% 95.5%
IC2RackLHS = [geomeanRackRatioLHS_/(multiStdRackLHS^2) geomeanRackRatioLHS_*(multiStdRackLHS^2)];
IC2RackRHS = [geomeanRackRatioRHS_/(multiStdRackRHS^2) geomeanRackRatioRHS_*(multiStdRackRHS^2)];
% 99.7%
IC3RackLHS = [geomeanRackRatioLHS_/(multiStdRackLHS^3) geomeanRackRatioLHS_*(multiStdRackLHS^3)];
IC3RackRHS = [geomeanRackRatioRHS_/(multiStdRackRHS^3) geomeanRackRatioRHS_*(multiStdRackRHS^3)];

% True IC representation
repIC1RackLHS = sum((omitZerosLHS < IC1RackLHS(2)) & (omitZerosLHS > IC1RackLHS(1)))/NLHS*100;
repIC1RackRHS = sum((omitZerosRHS < IC1RackRHS(2)) & (omitZerosRHS > IC1RackRHS(1)))/NRHS*100;
repIC2RackLHS = sum((omitZerosLHS < IC2RackLHS(2)) & (omitZerosLHS > IC2RackLHS(1)))/NLHS*100;
repIC2RackRHS = sum((omitZerosRHS < IC2RackRHS(2)) & (omitZerosRHS > IC2RackRHS(1)))/NRHS*100;
repIC3RackLHS = sum((omitZerosLHS < IC3RackLHS(2)) & (omitZerosLHS > IC3RackLHS(1)))/NLHS*100;
repIC3RackRHS = sum((omitZerosRHS < IC3RackRHS(2)) & (omitZerosRHS > IC3RackRHS(1)))/NRHS*100;

% Display summary statistics
fprintf('Rack demand ratio statistics (LHS, RHS):\n')
fprintf('Min: %.2f, %.2f\n',rackDemandStats(:,1))
fprintf('Max: %.2f, %.2f\n',rackDemandStats(:,2))
fprintf('Geometric mean: %.2f, %.2f\n',geomeanRackRatioLHS_,geomeanRackRatioRHS_)
fprintf('Multiplicative standard deviation: %.2f, %.2f\n',multiStdRackLHS,multiStdRackRHS)

% Display distribution limits
fprintf('Log-normal distribution limits of rack demand:\n')
fprintf('One standard deviation:\n')
disp([IC1RackLHS; IC1RackRHS])
fprintf('Two standard deviations:\n')
disp([IC2RackLHS; IC2RackRHS])
fprintf('Three standard deviations:\n')
disp([IC3RackLHS; IC3RackRHS])

% Display true confidence intervals
% fprintf('Log-normal rack demand confidence intervals:\n')
% fprintf('One standard deviation:\n')
% disp([repIC1RackLHS repIC1RackLHS])
% fprintf('Two standard deviations:\n')
% disp([repIC2RackLHS repIC2RackLHS])
% fprintf('Three standard deviations:\n')
% disp([repIC3RackLHS repIC3RackLHS])

repICRack = [repIC1RackLHS repIC1RackRHS;repIC2RackLHS repIC2RackRHS;repIC3RackLHS repIC3RackRHS];
fprintf('Confidence intervals using log-normal distribution properties:\n')
disp(repICRack)

% Portion of sample less than upper bound (captured by size recommendation)
lessThan1StdRack = [sum(omitZerosLHS < IC1RackLHS(2))/NLHS*100 sum(omitZerosRHS < IC1RackRHS(2))/NRHS*100];
lessThan2StdRack = [sum(omitZerosLHS < IC2RackLHS(2))/NLHS*100 sum(omitZerosRHS < IC2RackRHS(2))/NRHS*100];
lessThan3StdRack = [sum(omitZerosLHS < IC3RackLHS(2))/NLHS*100 sum(omitZerosRHS < IC3RackRHS(2))/NRHS*100];
lessThanStdRack = [lessThan1StdRack; lessThan2StdRack; lessThan3StdRack];

% Display cases captured
% fprintf('Rack demand cases captured by log-normal distribution properties:\n')
% fprintf('One standard deviation:\n')
% disp(lessThan1StdRack)
% fprintf('Two standard deviations:\n')
% disp(lessThan2StdRack)
% fprintf('Three standard deviations:\n')
% disp(lessThan3StdRack)

fprintf('Demand cases captured by log-normal distribution properties:\n')
disp(lessThanStdRack)

%% Rack demand ratio statistics
% Log-normal distribution (multiplicative normal distribution) by pitch

geomeanRackRatioArrayLHS = zeros(1,NPitchs);
geomeanRackRatioArrayRHS = zeros(1,NPitchs);

multiStdRackArrayLHS = zeros(1,NPitchs);
multiStdRackArrayRHS = zeros(1,NPitchs);

IC1RackLHS_ = zeros(2,NPitchs); 
IC1RackRHS_ = zeros(2,NPitchs);
IC2RackLHS_ = zeros(2,NPitchs);
IC2RackRHS_ = zeros(2,NPitchs);
IC3RackLHS_ = zeros(2,NPitchs);
IC3RackRHS_ = zeros(2,NPitchs);

for ii = 1:NPitchs
    
    andNowThisLHS = [];
    andNowThisRHS = [];
    
    for jj = 1:NARs
        
        andNowThisLHSInner = squeeze(rackDemandRatioLHSMatrix(:,jj,1:NTnsArray(jj),ii));
        andNowThisLHS = [andNowThisLHS andNowThisLHSInner(:)'];
        
        andNowThisRHSInner = squeeze(rackDemandRatioRHSMatrix(:,jj,1:NTnsArray(jj),ii));
        andNowThisRHS = [andNowThisRHS andNowThisRHSInner(:)']; 
        
    end
    
    % Omit zero results (required for log-normal distribution analysis)
    omitZerosLHS = andNowThisLHS(andNowThisLHS ~= 0);
    omitZerosRHS = andNowThisRHS(andNowThisRHS ~= 0);
    
    % Sample size
    NLHS = length(omitZerosLHS);
    NRHS = length(omitZerosRHS);
    
    % Geometric mean
    geomeanRackRatioArrayLHS(ii) = geomean(omitZerosLHS);
    geomeanRackRatioArrayRHS(ii) = geomean(omitZerosRHS);
    
    % Multiplicative standard deviation
    multiStdRackArrayLHS(ii) = exp(sqrt((sum((log(omitZerosLHS/geomeanRackRatioArrayLHS(ii))).^2))/(NLHS-1)));
    multiStdRackArrayRHS(ii) = exp(sqrt((sum((log(omitZerosRHS/geomeanRackRatioArrayRHS(ii))).^2))/(NRHS-1)));
    
    % Interval of confidence
    % 68.3%
    IC1RackLHS_(:,ii) = [geomeanRackRatioArrayLHS(ii)/multiStdRackArrayLHS(ii) geomeanRackRatioArrayLHS(ii)*multiStdRackArrayLHS(ii)];
    IC1RackRHS_(:,ii) = [geomeanRackRatioArrayRHS(ii)/multiStdRackArrayRHS(ii) geomeanRackRatioArrayLHS(ii)*multiStdRackArrayRHS(ii)];
    % 95.5%
    IC2RackLHS_(:,ii) = [geomeanRackRatioArrayLHS(ii)/(multiStdRackArrayLHS(ii)^2) geomeanRackRatioArrayLHS(ii)*(multiStdRackArrayLHS(ii)^2)];
    IC2RackRHS_(:,ii) = [geomeanRackRatioArrayRHS(ii)/(multiStdRackArrayRHS(ii)^2) geomeanRackRatioArrayLHS(ii)*(multiStdRackArrayRHS(ii)^2)];
    % 99.7%
    IC3RackLHS_(:,ii) = [geomeanRackRatioArrayLHS(ii)/(multiStdRackArrayLHS(ii)^3) geomeanRackRatioArrayLHS(ii)*(multiStdRackArrayLHS(ii)^3)];
    IC3RackRHS_(:,ii) = [geomeanRackRatioArrayRHS(ii)/(multiStdRackArrayRHS(ii)^3) geomeanRackRatioArrayLHS(ii)*(multiStdRackArrayRHS(ii)^3)];
    
end

geomeanRackRatioArray = [geomeanRackRatioArrayLHS; geomeanRackRatioArrayRHS];
multiStdRackArray = [multiStdRackArrayLHS; multiStdRackArrayRHS];
relativeGeomeanRackRatioArrayLHS = 100*geomeanRackRatioArrayLHS./geomeanRackRatioArrayLHS(1);
relativeGeomeanRackRatioArrayRHS = 100*geomeanRackRatioArrayRHS./geomeanRackRatioArrayRHS(1);
relativeGeomeanRackRatioArray = [relativeGeomeanRackRatioArrayLHS; relativeGeomeanRackRatioArrayRHS];

% Rack demand summary statistics by pitch
fprintf('Pitch size influence on rack demand ratio:\n')
fprintf('Geometric mean:\n')
disp(geomeanRackRatioArray)
fprintf('Multiplicative Standard deviation:\n')
disp(multiStdRackArray)
fprintf('Geometric mean as a percentage of 1mm pitch value:\n')
disp(relativeGeomeanRackRatioArray)

r2dDemandRatioLHSMatrix = rackDemandLHSMatrix./demandLHSMatrix;
r2dDemandRatioRHSMatrix = rackDemandRHSMatrix./demandRHSMatrix;

fprintf('Maximum rack demand to dissipater demand ratios:\n')
disp(max(r2dDemandRatioLHSMatrix(:)))
disp(max(r2dDemandRatioRHSMatrix(:)))

r2dAbove1LHS = length(find(r2dDemandRatioLHSMatrix > 1.0));
r2dAbove1RHS = length(find(r2dDemandRatioRHSMatrix > 1.0));

fprintf('Percentage of cases of greater demand in the rack than in the dissipater:\n')
disp(r2dAbove1LHS/Nanalyses*100)
disp(r2dAbove1RHS/Nanalyses*100)

return

%% Demand ratio (All)

for ii = 1:NARs
    
    legData = {};
    figure
    hold on
    for jj = 1:NTnsArray(ii)
        h = plot(Pitchs*1000,[squeeze(StatsArrayLHS(ii,jj,:,4)) squeeze(StatsArrayRHS(ii,jj,:,4))],'o--');
        set(h, {'marker'}, {'o'; 's'});
        myEntryLHS = strcat('Tn = ',num2str(TnsCell{ii}(jj)),' s',' (LHS)');
        myEntryRHS = strcat('Tn = ',num2str(TnsCell{ii}(jj)),' s',' (RHS)');
        legData = [legData myEntryLHS myEntryRHS];
    end
    grid on
    ylim([0 inf])
    xlabel('Pitch (mm)')
    ylabel('Geometric mean of demand ratio')
    legend(legData,'Location','Best')
    
end

%% Inelastic dissipater demand vs. peak uplift plots
% Break down by structure and pitch size

%array limits for demand ratio lines (from inspection)
rat1x = [0.4 0.4 0 0 0 0 0;...
        0.09 0.07 0.07 0.12 0 0 0;...
        0.18 0.16 0.15 0.18 0.15 0 0;...
        0.14 0.14 0.14 0.18 0.18 0.3 0.35];

rat10y = [1 1 0 0 0 0 0;...
        0.6 0.6 0.6 0.6 0 0 0;...
        0.7 0.7 0.7 0.7 0.7 0 0;...
        0.6 0.6 0.6 0.6 0.6 0.9 1.2];

% derived
rat1y = rat1x;
rat10x = rat10y/10;

for ii = 1:NARs
    
    for jj = 1:NTnsArray(ii)
        
        figure
        plot(squeeze(upliftLHS(:,ii,jj,1)),squeeze(demandLHSMatrix(:,ii,1)),'bo')
        grid on
        hold on
        plot(squeeze(upliftLHS(:,ii,jj,2)),squeeze(demandLHSMatrix(:,ii,jj,2)),'go')
        plot(squeeze(upliftLHS(:,ii,jj,3)),squeeze(demandLHSMatrix(:,ii,jj,3)),'ro')
        plot(squeeze(upliftLHS(:,ii,jj,4)),squeeze(demandLHSMatrix(:,ii,jj,4)),'co')
        plot(squeeze(upliftLHS(:,ii,jj,5)),squeeze(demandLHSMatrix(:,ii,jj,5)),'mo')
        plot([0 rat1x(ii,jj)],[0 rat1y(ii,jj)],'g--')
        plot([0 rat10x(ii,jj)],[0 rat10y(ii,jj)],'r--')
        legend('Pitch = 1mm','Pitch = 2mm','Pitch = 5mm','Pitch = 10mm',...
            'Pitch = 20mm','Ratio = 1','Ratio = 10','Location','Best')
        xlabel('Peak uplift (m)')
        ylabel('Inelastic dissipater demand (m)')
        title({'LHS dissipater demand ratios',...
            ['Period: ',num2str(TnsCell{ii}(jj)),' s, Aspect Ratio: ',num2str(ARs(ii))]})
        
        figure
        plot(squeeze(upliftRHS(:,ii,jj,1)),squeeze(demandRHSMatrix(:,ii,1)),'bo')
        grid on
        hold on
        plot(squeeze(upliftRHS(:,ii,jj,2)),squeeze(demandRHSMatrix(:,ii,jj,2)),'go')
        plot(squeeze(upliftRHS(:,ii,jj,3)),squeeze(demandRHSMatrix(:,ii,jj,3)),'ro')
        plot(squeeze(upliftRHS(:,ii,jj,4)),squeeze(demandRHSMatrix(:,ii,jj,4)),'co')
        plot(squeeze(upliftRHS(:,ii,jj,5)),squeeze(demandRHSMatrix(:,ii,jj,5)),'mo')
        plot([0 rat1x(ii,jj)],[0 rat1y(ii,jj)],'g--')
        plot([0 rat10x(ii,jj)],[0 rat10y(ii,jj)],'r--')
        legend('Pitch = 1mm','Pitch = 2mm','Pitch = 5mm','Pitch = 10mm',...
            'Pitch = 20mm','Ratio = 1','Ratio = 10','Location','Best')
        xlabel('Peak uplift (m)')
        ylabel('Inelastic dissipater demand (m)')
        title({'RHS dissipater demand ratios',...
            ['Period: ',num2str(TnsCell{ii}(jj)),' s, Aspect Ratio: ',num2str(ARs(ii))]})
        
    end
    
end

