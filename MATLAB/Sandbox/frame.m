function [ K, Khat, Lambda, KG ] = frame( E, A, I, L, theta, ne, ndof, Assem )
%FRAME Return stiffness properties for frame element
%   Detailed explanation goes here

%PREALLOCATION
K = zeros(6,6,ne);
Lambda = zeros(6,6,ne);
Khat = zeros(6,6,ne);
KG = zeros(ndof,ndof,ne);

for ii = 1:ne
    
    %LOCAL STIFFNESS
    alpha1 = E(ii)*I(ii)/(L(ii)^3);
    beta1 = (A(ii)*L(ii)^2)/I(ii);
    
    r1 = [beta1 0 0 -beta1 0 0];
    r2 = [0 12 6*L(ii) 0 -12 6*L(ii)];
    r3 = [0 6*L(ii) 4*(L(ii)^2) 0 -6*L(ii) 2*(L(ii)^2)];
    r4 = [-beta1 0 0 beta1 0 0];
    r5 = [0 -12 -6*L(ii) 0 12 -6*L(ii)];
    r6 = [0 6*L(ii) 2*(L(ii)^2) 0 -6*L(ii) 4*(L(ii)^2)];
    
    K(:,:,ii) = alpha1*[r1;r2;r3;r4;r5;r6];
    
    %GLOBAL STIFFNESS
    lambda = [cos(theta(ii)*pi/180) sin(theta(ii)*pi/180) 0;...
        -sin(theta(ii)*pi/180) cos(theta(ii)*pi/180) 0;0 0 1];
    Lambda(:,:,ii) = [lambda zeros(3);zeros(3) lambda];
    Khat(:,:,ii) = Lambda(:,:,ii)'*K(:,:,ii)*Lambda(:,:,ii);
    
    KG(:,:,ii) = Assem(:,:,ii)*Khat(:,:,ii)*Assem(:,:,ii)';
    
end

