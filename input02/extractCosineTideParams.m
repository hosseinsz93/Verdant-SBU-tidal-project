function [velMax, T, phi0] = extractCosineTideParams(DMY, velSigned)

    % Remove NaN values if need be
    if any(isnat(DMY)) && any(isnan(velSigned))
        velDMY_noNAN = ~isnat(DMY) & ~isnan(velSigned);
        DMY = DMY(velDMY_noNAN);
        velSigned = velSigned(velDMY_noNAN);
    end

    % Time converted to seconds from start
    t = seconds(DMY - DMY(1));

    % Tide's period is fixed as about semidiurnal (12 hrs & 25 mins according
    % to https://oceanservice.noaa.gov/facts/tidefrequency.html)
    T_fixed = 12*3600+25*60;

    % Guess made based on maxima in signed velocity
    velMax_guess = max(abs(velSigned));
    [~,indMax] = max(velSigned);
    phi0_guess = t(indMax);
    
    velMax_size = size(velMax_guess);
    phi0_size = size(phi0_guess);
    
    if velMax_size == phi0_size
        %fprintf("No transpose needed to correct array dimensions between maximum signed velocities and phi0.\n");
    elseif velMax_size(1) == phi0_size(2) && velMax_size(2) == phi0_size(1)
        %fprintf("Transpose needed to correct array dimensions between maximum signed velocities and phi0.\n");
        phi0_guess = phi0_guess.';
    else
        error("The array's length of maximum signed velocities (%i) is not equal to that of phi0 (%i).\n", length(velMax_guess),length(phi0_guess));
    end

    % Fitting to cosine by least squares
    paramsGuess = [velMax_guess, phi0_guess];
    cosineFit = @(params,t) params(1) * cos(2*pi(t-params(2))/T_fixed-pi/2);
    bounds = [0 0; 2*velMax_guess max(t)];

    opts = optimoptions('lsqcurvefit','Display','off');
    try
        paramsFit = lsqcurvefit(cosineFit, paramsGuess, t, velSigned, bounds(2,:), bounds(1,:), opts);
        velMax = paramsFit(1);
        phi0 = paramsFit(2);
        T = T_fixed;
    catch
        % If fitting fails, use initial guesses
        velMax = velMax_guess;
        phi0 = phi0_guess;
        T = T_fixed;
        warning('Cosine fitting failed, using initial parameter estimates');        
    end

end