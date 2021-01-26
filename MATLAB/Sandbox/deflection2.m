function [ defcoords ] = deflection2( coords, D, ne, te )
%DEFLECTION Produce coordinates of deformed node positions from original
%coordinates and FEA displacements

defcoords = zeros(ne,4,te);

for jj = 1:te
    
    defcoords(:,:,jj) = coords;
    
    for ii = 1:ne %element
        
        defcoords(ii,:,jj) = [coords(ii,1)+D(1,ii,jj) coords(ii,2)+D(2,ii,jj)...
            coords(ii,3)+D(4,ii,jj) coords(ii,4)+D(5,ii,jj)];
        
    end
    
end

end

