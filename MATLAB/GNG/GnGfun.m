function [ GnG ] = GnGfun( x, GnG, ii )
%GnG hysteresis behaviour model with ratcheting
%   Menegotto-Pinto material model with strain hardening
%   Returns corresponding input force for a given displacement and updates
%   reset and yield values for the algorithm

GnG.dx(ii) = x - GnG.xp; % DISPLACEMENT INPUT STEP
sx = sign(GnG.dx(ii)); % SIGN OF DISPLACEMENT INPUT STEP
GnG.fuse(ii) = GnG.fuse(ii-1); % UPDATE FUSE LENGTH (IN CASE OF NO CHANGE)

if GnG.fuse(ii) > GnG.fusebreak
    
    return
    
end

if x > GnG.x0 % IF BEYOND ENGAGEMENT DISPLACEMENT
    
%     if sx*GnG.sx_p < 0 % IF SIGN OF VELOCITY CHANGES
    
    if sx < 0 % IF SIGN OF VELOCITY CHANGES
        
        GnG.xr = GnG.xp; % UPDATE RESET DISPLACEMENT
        GnG.Fr = GnG.Fp; % UPDATE RESET FORCE
        GnG.x0 = GnG.xr - (GnG.Fr / GnG.K1); % UPDATE ENGAGEMENT DISPLACEMENT
        GnG.xy = GnG.x0 + GnG.Fy / GnG.K1; % UPDATE YIELD DISPLACEMENT
        
        %         if GnG.xp >= GnG.xy % IF QUASI-YIELDING HAS OCCURRED
        %
        %             GnG.xr = GnG.xp; % UPDATE RESET DISPLACEMENT
        %             GnG.Fr = GnG.Fp; % UPDATE RESET FORCE
        %             GnG.x0 = GnG.xr - (GnG.Fr / GnG.K1); % UPDATE ENGAGEMENT DISPLACEMENT
        %             GnG.xy = GnG.x0 + GnG.Fy / GnG.K1; % UPDATE YIELD DISPLACEMENT
        
        if GnG.Fp > GnG.Fy % IF YIELDING HAS OCCURRED
            
            GnG.Fy = GnG.Fp; % UPDATE YIELD FORCE
            GnG.xy = GnG.xp; % UPDATE YIELD DISPLACEMENT
            
        end
        
        %         end
        
    end
    
    % force (Menegotto-Pinto)
    GnG.F(ii) = (GnG.K1*(x - GnG.x0)) / ((1 + abs(GnG.K1*(x - GnG.x0) / GnG.Fy)^GnG.Beta)^(1/GnG.Beta));
    
    % original MP equation for reference
%     F = (K1*(x - xr) + Fr)/(1 + abs((K1*(x - xr) + Fr)/(Fy*sx))^Beta)^(1/Beta);
    
    if (x - GnG.xy) > 0 % if beyond yield
        
        GnG.F(ii) = GnG.F(ii) + GnG.K2*(x - GnG.xy); % allow for strain hardening (**EDIT: Fy for F**)
        
        if GnG.xp < GnG.xy % check if entire displacement step is yield
            
            GnG.fuse(ii) = GnG.fuse(ii-1) + GnG.dx(ii) - (GnG.xy - GnG.xp); % update fuse length
            
        else
            
            GnG.fuse(ii) = GnG.fuse(ii-1) + GnG.dx(ii); % update fuse length
            
        end
        
    end
    
    if GnG.dx(ii) < 0 % if unloading follow linear stiffness line
        
        GnG.F(ii) = GnG.K1*(x - GnG.xr) + GnG.Fr; % linear stiffness line
        
    end
    
else
    
    GnG.F(ii) = 0; % no force below engagement displacement
    
end

GnG.F(ii) = max([GnG.F(ii) 0]); % remove compressive forces

if length(GnG.ratchet) < GnG.teeth && x < GnG.x0 - GnG.pitch % check remaining...
    % ratchet capacity and compressive displacement
    
    GnG.ratchet = [GnG.ratchet, ii]; % record data point reference
    GnG.x0 = GnG.x0 - GnG.pitch; % adjust engagement displacement
    
    GnG.xy = GnG.x0 + (GnG.Fy / GnG.K1); % update yield displacement
    GnG.xr = GnG.x0 + (GnG.Fr / GnG.K1); % update reset displacement
    
end

GnG.sx_p = sx; % update sign of previous displacement step
GnG.xp = x; % update previous displacement
GnG.Fp = GnG.F(ii); % update previous force value

end
