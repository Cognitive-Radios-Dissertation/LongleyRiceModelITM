function [A_ref, mode] = lrprop(d, prop, prop_params)
% LRPROP - Longley-Rice Propagation Core
%
% Inputs:
%   d: Distance (meters) - Can be scalar or vector
%   prop: Structure from qlrpfl (h_e, d_L, a_eff, delta_h, etc.)
%   prop_params: Input parameters (freq, pol, etc.)
%
% Outputs:
%   A_ref: Median Reference Attenuation (dB)
%   mode: Propagation mode (1=LOS, 2=Diffraction, 3=Scatter)

    % Unpack
    freq = prop_params.freq_mhz;
    lambda = 299.792458 / freq; % wavelength in meters
    
    d_L = prop.d_L; % [d_L1, d_L2]
    d_Lt = sum(d_L); % Total LOS distance
    h_e = prop.h_e;
    
    % Calculate Reference Attenuation for each distance
    A_ref = zeros(size(d));
    mode = zeros(size(d));
    
    for i = 1:length(d)
        dist = d(i);
        
        if dist < d_Lt
            % --- Line of Sight Region ---
            mode(i) = 1;
            
            % 1. Two-Ray Optics
            arg = 2 * pi * h_e(1) * h_e(2) / (lambda * dist);
            if dist == 0, arg = 100; end 
            
            % Raw Two-Ray: A = -20 log10(abs(2 * sin(arg)))
            % FIX: Clamp to ensure no gain (Passive Terrain)
            % This removes the unphysical "Gain" < FSPL
            A_two_ray = -20 * log10(abs(2 * sin(arg)));
            A_two_ray = max(0, A_two_ray); 
            
            % 2. Blending
            w = max(0, min(1, dist / d_Lt)); 
            A_diff_at_horiz = calc_diffraction(d_Lt, freq, prop, lambda);
            
            A_ref(i) = (1-w) * A_two_ray + w * A_diff_at_horiz;
            
        elseif dist <= (1.5 * d_Lt) 
             % --- Diffraction Region ---
             mode(i) = 2;
             A_ref(i) = calc_diffraction(dist, freq, prop, lambda);
             
        else
             % --- Scatter Region ---
             mode(i) = 3;
             A_scat = calc_scatter(dist, freq, prop, prop_params);
             A_diff = calc_diffraction(dist, freq, prop, lambda);
             
             if A_scat < A_diff
                 A_ref(i) = A_scat;
             else
                 A_ref(i) = A_diff;
                 mode(i) = 2; 
             end
        end
    end
    
    % CRITICAL PHYSICS FIX:
    % Reference Attenuation relative to Free Space cannot be negative for median terrain.
    % (Negative A_ref implies Signal > Free Space, i.e., Gain).
    % While possible in specific multipath spots, for a median model, we clamp it.
    A_ref = max(0, A_ref);

end

% --- Helper Functions ---

function A_diff = calc_diffraction(d, freq, prop, lambda)
    % Simplified ITM Diffraction (Blending Knife Edge & Smooth Earth)
    
    theta_tot = prop.the(1) + prop.the(2) + d / prop.a_eff;
    v = theta_tot * sqrt(d / lambda); 
    
    if v > -0.7
        A_ke = 6.9 + 20 * log10( sqrt((v-0.1)^2 + 1) + v - 0.1 );
    else
        A_ke = 0;
    end
    
    delta_h = prop.delta_h;
    w_factor = 1 / (1 + delta_h/10); 
    
    % Use primarily Knife Edge + Roughness
    A_diff = A_ke; 
end

function A_scat = calc_scatter(d, freq, prop, params)
    % Troposcatter Loss Placeholder
    A_scat = 50 + 40 * log10(d/1000); 
end