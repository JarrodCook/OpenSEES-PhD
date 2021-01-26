function [ Assem ] = Assembly( dofmatch )
%ASSEMBLY Create assembly matrices from dof matches

Assem = zeros(max(max(dofmatch)),6,size(dofmatch,2));

dofref = eye(6);

for ii = 1:size(dofmatch,2) %element
    
    for jj = 1:6 %local DOF
        
        if dofmatch(jj,ii) ~= 0 %unfixed DOF
            
            Assem(dofmatch(jj,ii),:,ii)  = dofref(jj,:);
            
        end
        
    end
    
end

end

