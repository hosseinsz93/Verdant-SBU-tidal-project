function [floodVel, ebbVel, isFlood, isEbb, signedVel] = identifyTidalFlow(velocityMagnitude, direction)
    % IDENTIFYTIDALFLOW - Identifies flood and ebb tides from ADCP data
    %
    % This function classifies ADCP measurements into flood or ebb tide periods
    % based on the flow direction. The classification algorithm follows the 
    % approach used in TimeAveragedProfilesContinuous.m.
    %
    % Syntax:
    %   [floodVel, ebbVel, isFlood, isEbb, signedVel] = identifyTidalFlow(velocityMagnitude, direction)
    %
    % Inputs:
    %   velocityMagnitude - Vector of velocity magnitudes (m/s)
    %   direction        - Vector of flow directions (degrees, 0-360)
    %
    % Outputs:
    %   floodVel        - Velocity magnitudes during flood tide (m/s)
    %   ebbVel          - Velocity magnitudes during ebb tide (m/s)
    %   isFlood         - Boolean array indicating flood tide periods
    %   isEbb           - Boolean array indicating ebb tide periods
    %   signedVel       - Signed velocity (positive for flood, negative for ebb)
    %
    % Algorithm:
    %   - Directions 180°-360° = flood tide (eastward flow, positive sign)
    %   - Directions 0°-180° = ebb tide (westward flow, negative sign)
    %
    % Example:
    %   [floodVel, ebbVel, isFlood, isEbb, signedVel] = identifyTidalFlow(velocity, direction);
    %   floodAvg = mean(floodVel);
    %   ebbAvg = mean(ebbVel);
    %
    % Author: Hossein Seyedzadeh
    % Date: April 16, 2025
    
    % Input validation
    if length(velocityMagnitude) ~= length(direction)
        error('Velocity magnitude and direction vectors must be the same length');
    end
    
    % Initialize arrays
    n = length(velocityMagnitude);
    binaryDir = zeros(n, 1);
    signedVel = zeros(n, 1);
    isFlood = false(n, 1);
    isEbb = false(n, 1);
    
    % Step 1: Direction-based classification
    % Per TimeAveragedProfilesContinuous algorithm:
    % - 0° to 180° = ebb tide (westward flow, negative)
    % - 180° to 360° = flood tide (eastward flow, positive)
    for i = 1:n
        if (0 <= direction(i)) && (direction(i) < 180)
            binaryDir(i) = -1;  % West flow (ebb)
            isEbb(i) = true;
        elseif (180 <= direction(i)) && (direction(i) < 360)
            binaryDir(i) = 1;   % East flow (flood)
            isFlood(i) = true;
        end
    end
    
    % Step 2: Apply direction to velocity magnitude to get signed velocity
    % Positive velocity = flood tide, negative velocity = ebb tide
    signedVel = velocityMagnitude .* binaryDir;
    
    % Step 3: Extract velocity components
    % Create arrays containing only flood or ebb velocities
    floodVel = velocityMagnitude .* isFlood;  % Zero where not flood
    ebbVel = velocityMagnitude .* isEbb;      % Zero where not ebb
    
    % Remove zeros from the extracted velocities
    % This ensures statistics only include actual flood/ebb values
    floodVel = floodVel(isFlood);
    ebbVel = ebbVel(isEbb);
    
    % Optional: Identify major flood/ebb events using peak detection
    % Uncomment if you want this functionality
    %{
    velStd = std(signedVel, 'omitnan');
    [~, floodPeakLocs] = findpeaks(signedVel, 'MinPeakProminence', velStd);
    [~, ebbPeakLocs] = findpeaks(-signedVel, 'MinPeakProminence', velStd);
    
    % Create arrays identifying major events
    isFloodPeak = false(n, 1);
    isEbbPeak = false(n, 1);
    isFloodPeak(floodPeakLocs) = true;
    isEbbPeak(ebbPeakLocs) = true;
    %}
end