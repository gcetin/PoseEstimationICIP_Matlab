%% =========================================================
%% FULL ORIENTATION ERROR (RotErrMag) — MAIN METRIC
%% =========================================================

fprintf("\n=== FULL ORIENTATION ERROR (degrees) ===\n");

labels = ["Planar 2D Eight","Appended 2D Eight","3D Eight"];
nRuns  = numel(T);   % T is your cell array of tables

Mean   = zeros(nRuns,1);
Median = zeros(nRuns,1);
Std    = zeros(nRuns,1);
P95    = zeros(nRuns,1);
N      = zeros(nRuns,1);

for i = 1:nRuns
    Ti = T{i};

    % --- Ensure RotErrMag exists ---
    if ~any(strcmp(Ti.Properties.VariableNames,"RotErrMag"))
        dRoll  = wrapTo180(Ti.Est_Roll  - Ti.GT_Roll);
        dPitch = wrapTo180(Ti.Est_Pitch - Ti.GT_Pitch);
        dYaw   = wrapTo180(Ti.Est_Yaw   - Ti.GT_Yaw);
        Ti.RotErrMag = sqrt(dRoll.^2 + dPitch.^2 + dYaw.^2);
    end

    x = double(Ti.RotErrMag(:));
    x = x(~isnan(x));

    N(i)      = numel(x);
    Mean(i)   = mean(x);
    Median(i) = median(x);
    Std(i)    = std(x);
    P95(i)    = prctile(x,95);

    fprintf("%-20s | N=%3d | Mean=%6.2f | Median=%6.2f | Std=%6.2f | P95=%6.2f\n", ...
            labels(i), N(i), Mean(i), Median(i), Std(i), P95(i));
end

RotStats = table(labels(:), N, Mean, Median, Std, P95, ...
    'VariableNames', ...
    ["Trajectory","N","Mean_deg","Median_deg","Std_deg","P95_deg"]);

disp(RotStats);

writetable(RotStats, fullfile(outDir,"table_full_orientation_error_stats.csv"));
