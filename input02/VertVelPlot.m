%==========================================================================

% Goal: Plots of Current Station Info from The Race (LIS1001) from 7 to 9
% May, 2010
% Maker: Jonathan Craig
% Day: 10 February, 2025

%==========================================================================
%%
clc;
clear;
close all;

%% Import from Excel CSV Files

workfolder = pwd;
[files_num,VarKey,C,DMY,dateInterval,convertFactor,siteID,seadepths,depth_units] = FolderReadCSV(workfolder);
ncol = size(C,2);
nrow = size(C,1);
velMag = zeros(nrow,ncol);
dir = zeros(nrow,ncol);

% Alternative: Read from corresponding txt file
% E.g.: testread = fileread("C:\Users\...\Downloads\LIS1016.txt");

% Folder made to contain outputs
mainOutputDir = fullfile(workfolder, 'depth_averaged_results');
if ~exist(mainOutputDir, 'dir')
    mkdir(mainOutputDir);
    fprintf('Created main output folder: %s\n', mainOutputDir);
end

%% Preformatting

% VelStart and DirStart can be explicitly programmed as the formatting
% should be consistent; however, the code determines these indices from the
% variables to be safe

velDir_key = find(~contains(VarKey,"Date"));

for i = 1:numel(VarKey)
    if contains(VarKey(i),"Speed") && ~exist('VelStart','var')
        VelStart = find(velDir_key == i);
    elseif contains(VarKey(i),"Dir") && ~exist('DirStart','var')
        DirStart = find(velDir_key == i);
    end
end

velMag = C(:,VelStart:2:end);
dir = C(:,DirStart:2:end);

velMag = rmmissing(velMag,2);
dir = rmmissing(dir,2);
%% Data conversion

% Bathymetric data are given by NOAA Tides and Currents (https://tidesandcurrents.noaa.gov/cdata/StationInfo?id=LIS1001#)

for i = 1:length(seadepths)
    depth_disp = sprintf("Depth %d of %d: ",i,length(seadepths));
    % seadepths(i) = input(depth_disp); % For files without depth in
    % formatting
    if seadepths(i)>0
        seadepths(i) = -1*seadepths(i);
    end
end

seatop = max(seadepths);
seabottom = min(seadepths);
% depth_units = input("Please specify the depth units for plotting purposes: ",'s');

velMag = velMag.*convertFactor; % Conversion from knots to meters/second
velDir_depthAvg = mean(dir,2,'omitnan');

% Calculation of velocity by introducing direction into magnitude
Eas = velMag.*cosd(dir);
Eas_depthAvg = mean(Eas,2,'omitnan');

Nor = velMag.*sind(dir);
Nor_depthAvg = mean(Nor,2,'omitnan');

%W = zeros(size(Eas));

%% Depth averaging

[floodVelMag, ebbVelMag, floodVelDir, ebbVelDir, floodVelMag_depthAvg, floodVelDir_depthAvg, ...
    ebbVelMag_depthAvg, ebbVelDir_depthAvg, floodSign, ebbSign, velMagSigned, ebbSign_depthAvg, floodSign_depthAvg] = identifyTidalFlow(velMag,dir);

%% Principal flow direction

floodVelDir_depthAvg_Trig = sin(2*floodVelDir_depthAvg)./cos(2*floodVelDir_depthAvg);
ebbVelDir_depthAvg_Trig = sin(2*ebbVelDir_depthAvg)./cos(2*ebbVelDir_depthAvg);

floodVelDir_depthtimeAvg_Trig = mean(floodVelDir_depthAvg_Trig,'omitnan');
ebbVelDir_depthtimeAvg_Trig = mean(ebbVelDir_depthAvg_Trig,'omitnan');

floodMeanDir = 0.5*atand(floodVelDir_depthtimeAvg_Trig);
ebbMeanDir = 0.5*atand(ebbVelDir_depthtimeAvg_Trig);

% Time-average of depth-averaged flood and ebb magnitude taken for flow
% visualization
floodVelMag_depthtimeAvg = mean(floodVelMag_depthAvg,'omitnan'); 
ebbVelMag_depthtimeAvg = mean(ebbVelMag_depthAvg,'omitnan');

% Alternative based on ebbFloodClassifier.ipynb

AB = horzcat(Eas_depthAvg,Nor_depthAvg);
AB_centered = AB-mean(AB,1);
AB_cov = cov(AB(:,1),AB(:,end));
[AB_eigvec,AB_eigval] = eig(AB_cov);
[AB_eigvalMax,imax] = max(AB_eigval);
eigVec = AB_eigvec(:,imax);
phi_deg = atan2d(eigVec(1),eigVec(2));
% mod(atan2d(,),360);

velStream = sind(phi_deg) + cosd(phi_deg);
velSpan = cosd(phi_deg) - sind(phi_deg);
%velPhase = 0;


%% Cosine Tidal Params

velSigned_depthAvg = mean(velMagSigned,2); % Converted to m/s
[velMax, T, phi0] = extractCosineTideParams(DMY, velSigned_depthAvg);
    
% Printed info about fitting cosine
fprintf("Cosine fit parameters:\n U_max = %.3f m/s\n T = %.1f s (%.2f h)\n phi0 = %.1f s (%.2f h)\n", ... 
    velMax, T, T/3600, phi0, phi0/3600);
fprintf("\n Fitted Cosine for Depth-Averaged Velocity:\n");
fprintf('    u(t) = %.3f * cos(2\pi*(t - %.1f)/%.1f - \pi/2)\n', velMax, phi0, T);
fprintf('    where t is time in seconds from start of record\n');
fprintf('    Alternative form: u(t) = %.3f * sin(2π*(t - %.1f)/%.1f)\n', velMax, phi0, T);

% Building the function and the fitting stats
t = seconds(DMY - DMY(1));
velCosFit = velMax*cos(2*pi*(t-phi0)/T-pi/2);
residuals = velMagSigned - velCosFit;
rmse = sqrt(mean(residuals.^2,'omitnan'));
r_squared = 1 - sum(residuals.^2, 'omitnan') / sum((velMagSigned - mean(velMagSigned, 'omitnan')).^2, 'omitnan');
fprintf('Fit quality: RMSE = %.3f m/s, R² = %.3f\n', rmse, r_squared);

% Plot of cosine fit
% Plot observed vs. fitted cosine
figure('Name', sprintf('%s (%s to %s): Cosine Fit',siteID,DMY(1),DMY(end)), 'Position', [100, 100, 1200, 400]);
plot(DMY, velSigned_depthAvg, 'k-', 'DisplayName', 'Observed');
hold on;
plot(DMY, velCosFit, 'r--', 'DisplayName', 'Cosine Fit');
plot(DMY, residuals, 'g:', 'DisplayName', 'Residuals');

legend;
xlabel('Time');
ylabel('Velocity (cm/s)');
title(sprintf('%s: Observed vs. Cosine-Fit Tidal Velocity from %s to %s',siteID,DMY(1),DMY(end)));
grid on;

% Add text box with fit statistics
text_str = sprintf('RMSE: %.3f m/s\nR²: %.3f\nMax Error: %.3f m/s', ...
    rmse, r_squared, max(abs(residuals), [], 'omitnan'));
text(0.02, 0.98, text_str, 'Units', 'normalized', 'VerticalAlignment', 'top', ...
    'BackgroundColor', 'white', 'EdgeColor', 'black', 'FontSize', 10);

saveas(gcf, fullfile(mainOutputDir, '_CosineFit.png'));

%% Turbulence Intensity

% Alternative data-filtering for ebb, flood, AND slack



%% File Storage

flowVisualized(mainOutputDir,DMY,velSigned_depthAvg,velDir_depthAvg,floodSign_depthAvg,ebbSign_depthAvg,Eas_depthAvg,Nor_depthAvg,floodVelMag_depthtimeAvg,ebbVelMag_depthtimeAvg,floodMeanDir,ebbMeanDir,siteID);

% flowSaved();

%% Profile characteristics

% siteProfile_size = [1 6];
% siteProfile_varTypes = ["string" "string" "double" "double" "double" "double"];
% siteProfile_varNames = ["Site code" "Tidal cycle type" "Misalignment between flood and ebb [/it{deg}]" "Peak speed [/it{m/s} (flood/ebb)]" "% duration greater than 4 m/s (flood/ebb)" "Reynolds number range [x 10^{7}] (LB-flood/ebb)"];
% siteProfile = table(Size=siteProfile_size, VariableTypes=siteProfile_varTypes, VariableNames=siteProfile_varNames);

