function [ xxdd, yydd ] = eulersf3( xx, yy, d, L, angle, shaperes, te )
%EULERSF Returns element displacement between nodes using Euler beam shape
%functions
%   Divides element into number of points defined by shaperes and
%   calculates the flexural and axial displacement throughout the element
%   via Euler beam shape functions

xxdd = zeros(shaperes,size(xx,2),te);
yydd = zeros(shaperes,size(xx,2),te);
dispfl = zeros(shaperes,size(xx,2),te);
dispax = zeros(shaperes,size(xx,2),te);

for jj = 1:te
    
    for ii = 1:size(xx,2) %element
        
        %interpolate element initial position
        xxdd(:,ii,jj) = linspace(xx(1,ii),xx(2,ii),shaperes);
        yydd(:,ii,jj) = linspace(yy(1,ii),yy(2,ii),shaperes);
        x = 0:L(ii)/(shaperes-1):L(ii);
        
        %FLEXURAL
        psi2 = 1-3/(L(ii)^2)*x.^2 + 2/(L(ii)^3)*(x.^3); %shape functions
        psi3 = x.*((1-x./L(ii)).^2);
        psi5 = 3/(L(ii)^2)*(x.^2) - 2/(L(ii)^3)*(x.^3);
        psi6 = (x.^3)/(L(ii)^2) - (x.^2)/L(ii);
        
        dispfl(:,ii,jj) = d(2,ii,jj)*psi2 + d(3,ii,jj)*psi3 + d(5,ii,jj)*psi5 + d(6,ii,jj)*psi6;
        %incorporate shape into interpolated vector
        xxdd(:,ii,jj) = xxdd(:,ii,jj) - dispfl(:,ii,jj)*sind(angle(ii));
        yydd(:,ii,jj) = yydd(:,ii,jj) + dispfl(:,ii,jj)*cosd(angle(ii));
        
        %AXIAL
        psi1 = 1 - x./L(ii); %shape functions
        psi2 = x./L(ii);
        
        dispax(:,ii,jj) = d(1,ii,jj)*psi1 + d(4,ii,jj)*psi2;
        %incorporate shape into interpolated vector
        xxdd(:,ii,jj) = xxdd(:,ii,jj) + dispax(:,ii,jj)*cosd(angle(ii));
        yydd(:,ii,jj) = yydd(:,ii,jj) + dispax(:,ii,jj)*sind(angle(ii));
        
    end
    
end

end

