function [ F, f, D, d ] = fordis2( K, Khat, Lambda, q, Assem, ne, te )
%FORDIS Returns forces and displacements in global and local coordinates
%for each element

F = zeros(6,ne,te);
D = zeros(6,ne,te);
d = zeros(6,ne,te);
f = zeros(6,ne,te);

for jj = 1:te
    
    for ii = 1:ne
        
        F(:,ii,jj) = Khat(:,:,ii)*Assem(:,:,ii)'*q(:,jj);
        D(:,ii,jj) = Assem(:,:,ii)'*q(:,jj);
        d(:,ii,jj) = Lambda(:,:,ii)*Assem(:,:,ii)'*q(:,jj);
        f(:,ii,jj) = K(:,:,ii)*d(:,ii,jj);
        
    end
    
end

end

