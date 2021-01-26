% Generic Linear Dynamic Setup
% jco118

close all; clear; clc; %reset
set(0,'units','centimeters','defaultfigureposition',[60 60 2*420,2*330]); %new default figure
tic %start the clock

%% COORDINATES

% Enter node coordinates for each element
% FORMAT: 1 row per element - [Xa Ya Xb Yb]

ne = 15; % # elements
coords = [0 0 0 1;0 1 1 1;1 0 1 1;0 1 0 2;0 2 1 2;1 1 1 2;1 2 2 2;2 1 2 2;1 1 2 1;...
    1 2 1 3;1 3 2 3;2 2 2 3;1 3 1 4;1 4 2 4;2 3 2 4];

xx = [coords(:,1)';coords(:,3)'];
yy = [coords(:,2)';coords(:,4)'];

%% ASSEMBLY

% Match element degrees of freedom (in global coordinates) to global
% degrees of freedom. 0 indicates fixed degree of freedom.

dofmatch = zeros(6,ne);

dofmatch(:,1) = [0 0 0 1 2 3];
dofmatch(:,2) = [1 2 3 4 5 6];
dofmatch(:,3) = [0 0 0 4 5 6];
dofmatch(:,4) = [1 2 3 7 8 9];
dofmatch(:,5) = [7 8 9 10 11 12];
dofmatch(:,6) = [4 5 6 10 11 12];
dofmatch(:,7) = [10 11 12 13 14 15];
dofmatch(:,8) = [16 17 18 13 14 15];
dofmatch(:,9) = [4 5 6 16 17 18];
dofmatch(:,10) = [10 11 12 19 20 21];
dofmatch(:,11) = [19 20 21 22 23 24];
dofmatch(:,12) = [13 14 15 22 23 24];
dofmatch(:,13) = [19 20 21 25 26 27];
dofmatch(:,14) = [25 26 27 28 29 30];
dofmatch(:,15) = [22 23 24 28 29 30];

[ Assem ] = Assembly( dofmatch );
ndof = max(dofmatch(:));

%% PROPERTIES

E = 200e9*ones(ne,1);
A = 1e-4*ones(ne,1);
I = 4e-6*ones(ne,1);
[ theta, L ] = anglefun( coords );

loadap = zeros(ndof,1);
loadap(25) = 1;
loadap(28) = 1;

%% STIFFNESS

[ K, Khat, Lambda, KG ] = frame( E, A, I, L, theta, ne, ndof, Assem );
KGsum = sum(KG,3);

%% ROTATIONAL SPRING

% Ktheta = 2e5;
% Asrs = [0 0; 0 0; 1 0; 0 1; 0 0; 0 0; 0 0];
% Khatrs = Ktheta*[1 -1; -1 1];
% KGrs = Asrs*Khatrs*Asrs';
% KGsum = KGsum + KGrs;

%% DYNAMIC PROPERTIES

M = diag(ones(ndof,1));
M(25,25) = 2e4;
M(28,28) = 2e4;
resfreq = sqrt(eig(KGsum,M));
coeffs = 0.5*[1/resfreq(1) resfreq(1);1/resfreq(2) resfreq(2)]\[0.05;0.05];
C = real(coeffs(1)*M + coeffs(2)*KGsum);

%% TIME AND FORCE

et = 30.0; % END TIME OF SIGNAL
dt = 0.04; % TIME STEP
te = et/dt + 1; % NUMBER OF TIME STEPS
t = linspace(0,et,te); % TIME ARRAY

acc = linspace(1,12.2,te).*(0.2*sin(1.5:1.5:1.5*te) + sin(0.2:0.2:0.2*te));
% acc = sin(t); % FORCING FUNCTION -- acceleration
dacc = [0 diff(acc)]; % FORCE INCREMENT ARRAY -- acceleration
p0 = M*loadap*acc; % INPUT LOAD FORCE
dp = M*loadap*dacc; % INPUT LOAD FORCE INCREMENT

%% NUMERICAL INTEGRATION

%PREALLOCATION
q = zeros(ndof,et);
xd = zeros(ndof,et);
xdd = zeros(ndof,et);

for ii = 2:te
    
    % INCREMENTAL NEWMARK BETA INTEGRATION
    [q(:,ii), xd(:,ii), xdd(:,ii-1)] = dynamic1(p0(:,ii-1), ...
        dp(:,ii), C, KGsum, q(:,ii-1), xd(:,ii-1), M, dt);
    
    e(:,ii) = M*xdd(:,ii-1) + C*xd(:,ii-1) + KGsum*q(:,ii-1) - p0(:,ii-1);
    
end

toc; %time

%% FORCES AND DISPLACEMENTS

[ F, f, D, d ] = fordis2( K, Khat, Lambda, q, Assem, ne, te );

[ defcoords ] = deflection2( coords, D, ne, te );

xxd = [permute(defcoords(:,1,:),[2 1 3]);permute(defcoords(:,3,:),[2 1 3])];
yyd = [permute(defcoords(:,2,:),[2 1 3]);permute(defcoords(:,4,:),[2 1 3])];

% xxd = [defcoords(:,1,:)';defcoords(:,3,:)'];
% yyd = [defcoords(:,2,:)';defcoords(:,4,:)'];

%% SHAPE FUNCTIONS

shaperes = 101; %shape function resolution
[ xxdd, yydd ] = eulersf3( xx, yy, d, L, theta, shaperes, te );

%% PLOTS

% figure
% h1 = plot(xx,yy,'k.');
% hold on;
% h2 = plot(xx,yy,'b-');
% h3 = plot(xxd,yyd,'r.');
% h4 = plot(xxdd,yydd,'g-');
% axis([min(xx(:))-1 max(xx(:))+1 min(yy(:))-1 max(yy(:))+1])
% axis equal
% legend([h1(1) h2(1) h3(1) h4(1)],'Undeflected (node)','Undeflected (element)',...
%     'Deflected (node)','Deflected (element)','Location','Best');
% 
% display(xxd)
% display(xxdd([1 shaperes],:))
% display(yyd)
% display(yydd([1 shaperes],:))

% Mov = zeros(te,1);

for ii = 1:te
    
    h = figure(1);
    
    plot(xxd(:,:,ii),yyd(:,:,ii),'r.',xxdd(:,:,ii),yydd(:,:,ii),'g-');
    axis([min(xx(:))-1 max(xx(:))+1 min(yy(:))-1 max(yy(:))+1])
    title('Structure Displacement')
    
%     Mov(ii) = getframe;
    
end

% movie(Mov,1);

% for ii = 1:te
%     
%     h = figure(1);
%     h1 = plot(xx,yy,'k.');
%     hold on;
%     h2 = plot(xx,yy,'b-');
%     h3 = plot(xxd(:,:,ii),yyd(:,:,ii),'r.');
%     h4 = plot(xxdd(:,:,ii),yydd(:,:,ii),'g-');
%     axis([min(xx(:))-1 max(xx(:))+1 min(yy(:))-1 max(yy(:))+1])
%     legend([h1(1) h2(1) h3(1) h4(1)],'Undeflected (node)','Undeflected (element)',...
%         'Deflected (node)','Deflected (element)','Location','Best');
%     
% end