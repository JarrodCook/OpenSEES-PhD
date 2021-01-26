function [ angle, L ] = anglefun( coords )
%ANGLEFUN Returns element angles and lengths from coordinates
% Enter node coordinates for each element
% FORMAT: 1 row per element - [Xa Ya Xb Yb]

angle = zeros(size(coords,1),1);
L = zeros(size(coords,1),1);

for ii = 1:size(coords,1)
    
    dx = coords(ii,3) - coords(ii,1);
    dy = coords(ii,4) - coords(ii,2);
    L(ii) = sqrt(dx^2+dy^2);
    
    if dx == 0
        angle(ii) = sign(dy)*90;
    elseif dy == 0
        angle(ii) = 90*(sign(dx)-1);
    elseif dx > 0
        angle(ii) = atand(dy/dx);
    else
        angle(ii) = atand(dy/dx) + 180*sign(dy);
    end
    
end

end

