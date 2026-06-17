function [floodVelMag, ebbVelMag, floodVelDir, ebbVelDir, floodVelMag_depthAvg, floodVelDir_depthAvg, ebbVelMag_depthAvg, ebbVelDir_depthAvg, floodSign, ebbSign, velMagSigned, ebbSign_depthAvg, floodSign_depthAvg] = identifyTidalFlow(velMag,dir)
% 
% Algorithm:
%  - Directions 180°-360° = flood tide (eastward flow, positive sign)
%  - Directions 0°-180° = ebb tide (westward flow, negative sign)
% 
% Original Author: Hossein Seyedzadeh
% Date Modified: 16 May, 2026

    % Input validation
    if size(velMag) ~= size(dir)
        error('Velocity magnitude and direction vectors must be the same size');
    end
    
    % Initialized arrays
    [rowMag,colMag] = size(velMag);
    binaryDir = zeros(rowMag,colMag);
    velMagSigned = zeros(rowMag,colMag);
    floodSign = false(rowMag,colMag);
    ebbSign = false(rowMag,colMag);
    
    for i = 1:rowMag
        for j = 1:colMag
            if (0<=dir(i,j)) && (dir(i,j)<180) || (dir(i,j)==360) % Inclusive of 360 degrees which represent full rotation to zero
                binaryDir(i,j) = -1;
                ebbSign(i,j) = true;
            elseif (180<=dir(i,j)) && (dir(i,j)<360)
                binaryDir(i,j) = 1;
                floodSign(i,j) = true;
            end
        end
    end
    
    % Exclusive vectors for ebb and flood
    floodVelMag = velMag.*floodSign;
    ebbVelMag = velMag.*ebbSign;
    
    floodVelDir = dir.*floodSign;
    ebbVelDir = dir.*ebbSign;
    
    % Preparation for principal flow
    floodVelMag_depthAvg = zeros(size(floodSign,1),1);
    floodVelDir_depthAvg = zeros(size(floodSign,1),1);
    ebbVelMag_depthAvg = zeros(size(ebbSign,1),1);
    ebbVelDir_depthAvg = zeros(size(ebbSign,1),1);
    floodSign_depthAvg = zeros(size(ebbSign,1),1);
    ebbSign_depthAvg = zeros(size(ebbSign,1),1);
    floodVelCount = 0;
    ebbVelCount = 0;

    for i = 1:rowMag
        for j = 1:colMag
            if floodSign(i,j)
                floodVelMag_depthAvg(i) = floodVelMag_depthAvg(i)+velMag(i,j);
                floodVelDir_depthAvg(i) = floodVelDir_depthAvg(i)+dir(i,j);
                floodVelCount = floodVelCount+1;
            elseif ebbSign(i,j)
                ebbVelMag_depthAvg(i) = ebbVelMag_depthAvg(i)+velMag(i,j);
                ebbVelDir_depthAvg(i) = ebbVelDir_depthAvg(i)+dir(i,j);
                ebbVelCount = ebbVelCount+1;
            else
                warning("Tide at i = %i and j = %i is marked as neither ebb nor flood.\n",i,j);
                continue
            end
        end
        
        floodVelMag_depthAvg(i) = floodVelMag_depthAvg(i)/floodVelCount;
        floodVelDir_depthAvg(i) = floodVelDir_depthAvg(i)/floodVelCount;
        floodSign_depthAvg(i) = ~isnan(floodVelDir_depthAvg(i));
        
        ebbVelMag_depthAvg(i) = ebbVelMag_depthAvg(i)/ebbVelCount;
        ebbVelDir_depthAvg(i) = ebbVelDir_depthAvg(i)/ebbVelCount;
        ebbSign_depthAvg(i) = ~isnan(ebbVelDir_depthAvg(i));

        % Reset count for next row
        floodVelCount = 0;
        ebbVelCount = 0;

    end
    
    floodSign_depthAvg = logical(floodSign_depthAvg);
    ebbSign_depthAvg = logical(ebbSign_depthAvg);
    
    % Preparation for fit to cosine
    velMagSigned = velMag.*binaryDir;
    
    % Remove zeros
%     floodVelMag = floodVelMag(floodSign);
%     ebbVelMag = ebbVelMag(ebbSign);

    % Optional: Identify major flood/ebb events using peak detection
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