clear; clc; close all;

files = {
".\Logs\DataFolder_1920x1080_50_NOISY_REL_MOTION_TRANS_XYZ_PAN_3D_EIGHT2D_5000it_1e-06minCost_250maxTrial_ang5_trans2_2_5.txt"
".\Logs\DataFolder_1920x1080_50_NOISY_REL_MOTION_TRANS_XYZ_PAN_3D_APP_EIGHT2D_5000it_1e-06minCost_250maxTrial_ang5_trans2_2_5.txt"
".\Logs\DataFolder_1920x1080_50_NOISY_REL_MOTION_TRANS_XYZ_PAN_3D_EIGHT3D_5000it_1e-06minCost_250maxTrial_ang5_trans2_2_5.txt"
};
labels = ["2D Eight","Appended 2D Eight","3D Eight"];

outDir = "analysis_outputs_compare";
if ~exist(outDir,"dir"); mkdir(outDir); end

% ----- Load tables
Ts = cell(numel(files),1);
for i=1:numel(files)
    opts = detectImportOptions(files{i},"Delimiter","\t");
    opts = setvartype(opts,"double");
    Ti = readtable(files{i},opts);

    if any(strcmp(Ti.Properties.VariableNames,"Success"))
        Ti = Ti(Ti.Success==1,:);
    end

    % Ensure RotErrMag exists
    if ~any(strcmp(Ti.Properties.VariableNames,"RotErrMag"))
        dRoll  = wrapTo180(Ti.Est_Roll  - Ti.GT_Roll);
        dPitch = wrapTo180(Ti.Est_Pitch - Ti.GT_Pitch);
        dYaw   = wrapTo180(Ti.Est_Yaw   - Ti.GT_Yaw);
        Ti.RotErrMag = sqrt(dRoll.^2 + dPitch.^2 + dYaw.^2);
    end

    % Translation error magnitude
    Ti.tErrMag = sqrt((Ti.Est_x-Ti.GT_x).^2 + (Ti.Est_y-Ti.GT_y).^2 + (Ti.Est_z-Ti.GT_z).^2);

    Ts{i} = Ti;
end

% ----- Choose binning variable (example: GT_z)
numBins = 12; % fewer bins for clearer grouped bars

% Common edges across all runs
allGTz = vertcat(Ts{1}.GT_z, Ts{2}.GT_z, Ts{3}.GT_z);
edges   = linspace(min(allGTz), max(allGTz)+eps, numBins+1);
centers = 0.5*(edges(1:end-1)+edges(2:end));

% ----- Compute mean per bin for each trajectory
M_rot = nan(numBins, numel(Ts));
M_trn = nan(numBins, numel(Ts));

for i=1:numel(Ts)
    Ti = Ts{i};
    b  = discretize(Ti.GT_z, edges);

    valid = ~isnan(b);
    M_rot(:,i) = accumarray(b(valid), Ti.RotErrMag(valid), [numBins 1], @mean, NaN);
    M_trn(:,i) = accumarray(b(valid), Ti.tErrMag(valid),   [numBins 1], @mean, NaN);
end

% ----- Plot: RotErrMag vs GT_z (grouped)
figure('Color','w');
bar(centers, M_rot, 1.0); grid on;
xlabel("GT Z (binned)");
ylabel("Mean RotErrMag (deg)");
title("Rotation error vs GT Z (binned) - trajectory comparison");
legend(labels, "Location","northwest");

% Sparse ticks
nt = min(7, numBins);
tickIdx = round(linspace(1,numBins,nt));
xticks(centers(tickIdx));
xticklabels(string(round(centers(tickIdx),1)));

saveas(gcf, fullfile(outDir,"RotErr_vs_GTz_grouped.png"));

% ----- Plot: tErrMag vs GT_z (grouped)
figure('Color','w');
bar(centers, M_trn, 1.0); grid on;
xlabel("GT Z (binned)");
ylabel("Mean translation error magnitude");
title("Translation error vs GT Z (binned) - trajectory comparison");
legend(labels, "Location","northwest");

xticks(centers(tickIdx));
xticklabels(string(round(centers(tickIdx),1)));

saveas(gcf, fullfile(outDir,"tErrMag_vs_GTz_grouped.png"));
