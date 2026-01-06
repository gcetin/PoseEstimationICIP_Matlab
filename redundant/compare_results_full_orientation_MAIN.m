clear; clc; close all;

%% ------------------ FILES ------------------
files = {
    ".\Logs\DataFolder_1920x1080_50_NOISY_REL_MOTION_TRANS_XYZ_PAN_3D_EIGHT2D_5000it_1e-06minCost_250maxTrial_ang5_trans2_2_5.txt"      % Planar 2D Eight
    ".\Logs\DataFolder_1920x1080_50_NOISY_REL_MOTION_TRANS_XYZ_PAN_3D_APP_EIGHT2D_5000it_1e-06minCost_250maxTrial_ang5_trans2_2_5.txt"  % Appended 2D Eight
    ".\Logs\DataFolder_1920x1080_50_NOISY_REL_MOTION_TRANS_XYZ_PAN_3D_EIGHT3D_5000it_1e-06minCost_250maxTrial_ang5_trans2_2_5.txt"      % 3D Eight
};
labels = ["Planar 2D Eight","Appended 2D Eight","3D Eight"];

outDir = "analysis_outputs_compare_fullOri";
if ~exist(outDir,"dir"); mkdir(outDir); end

numBins = 12;                 % for GT-z binned plots
binBy   = "GT_z";             % recommended
saveFigs = true;

%% ------------------ LOAD + METRICS ------------------
T = cell(numel(files),1);
for i=1:numel(files)
    T{i} = readtable(files{i}, detectImportOptions(files{i}, "Delimiter","\t"));
    if any(strcmp(T{i}.Properties.VariableNames,"Success"))
        T{i} = T{i}(T{i}.Success==1,:);
    end

    % Translation magnitude error
    T{i}.tErrMag = sqrt((T{i}.Est_x - T{i}.GT_x).^2 + (T{i}.Est_y - T{i}.GT_y).^2 + (T{i}.Est_z - T{i}.GT_z).^2);

    % FULL orientation magnitude (MAIN) — use RotErrMag if present, else compute proxy
    if ~any(strcmp(T{i}.Properties.VariableNames,"RotErrMag"))
        dR = wrapTo180(T{i}.Est_Roll  - T{i}.GT_Roll);
        dP = wrapTo180(T{i}.Est_Pitch - T{i}.GT_Pitch);
        dY = wrapTo180(T{i}.Est_Yaw   - T{i}.GT_Yaw);
        T{i}.RotErrMag = sqrt(dR.^2 + dP.^2 + dY.^2);
        warning("RotErrMag missing in %s → Euler L2 proxy computed.", labels(i));
    end

    fprintf("[Loaded] %-18s | N=%d\n", labels(i), height(T{i}));
end

%% ------------------ PRINT STATS (MAIN RESULTS) ------------------
RotStats = statsTable(T, labels, "RotErrMag");
TransStats = statsTable(T, labels, "tErrMag");

fprintf("\n=== FULL ORIENTATION ERROR (RotErrMag, degrees) ===\n");
disp(RotStats);

fprintf("\n=== TRANSLATION ERROR MAGNITUDE (tErrMag, GT units e.g., meters) ===\n");
disp(TransStats);

writetable(RotStats,  fullfile(outDir,"table_full_orientation_stats.csv"));
writetable(TransStats,fullfile(outDir,"table_translation_stats.csv"));

%% ------------------ FIGURE: FULL ORIENTATION vs GT_z (BINNED, GROUPED BAR) ------------------
[edges, centers] = commonBins(T, binBy, numBins);
[Ym, Ys] = binnedMeanStd(T, binBy, edges, "RotErrMag");

figure('Color','w'); grid on; hold on;
bar(centers, Ym, 1.0);
xlabel("GT Z (binned)");
ylabel("Mean orientation error magnitude (degrees)");
title("Full orientation error vs GT Z (binned)");
legend(labels,'Location','northwest');

% error bars = std (optional)
bh = flipud(findobj(gca,'Type','Bar'));
for k=1:numel(bh)
    errorbar(bh(k).XEndPoints, bh(k).YEndPoints, Ys(:,k), 'k.', 'LineWidth',1, 'CapSize',3);
end
hold off;

if saveFigs
    saveas(gcf, fullfile(outDir,"fig_fullOri_vs_GTz_binned_grouped.png"));
end

%% ------------------ FIGURE: TRANSLATION vs GT_z (BINNED, GROUPED BAR) ------------------
[YmT, YsT] = binnedMeanStd(T, binBy, edges, "tErrMag");

figure('Color','w'); grid on; hold on;
bar(centers, YmT, 1.0);
xlabel("GT Z (binned)");
ylabel("Mean translation error magnitude (GT units)");
title("Translation error vs GT Z (binned)");
legend(labels,'Location','northwest');

bh = flipud(findobj(gca,'Type','Bar'));
for k=1:numel(bh)
    errorbar(bh(k).XEndPoints, bh(k).YEndPoints, YsT(:,k), 'k.', 'LineWidth',1, 'CapSize',3);
end
hold off;

if saveFigs
    saveas(gcf, fullfile(outDir,"fig_trans_vs_GTz_binned_grouped.png"));
end

%% ------------------ FIGURE: CDFs (FULL ORIENTATION + TRANSLATION) ------------------
plotCDF(T, labels, "RotErrMag", "Orientation error magnitude (degrees)", ...
        "Full orientation error distribution (CDF)", fullfile(outDir,"fig_fullOri_cdf.png"), saveFigs);

plotCDF(T, labels, "tErrMag", "Translation error magnitude (GT units)", ...
        "Translation error distribution (CDF)", fullfile(outDir,"fig_trans_cdf.png"), saveFigs);

fprintf("\n✔ Done. Outputs saved to: %s\n", outDir);

%% ====================== HELPERS ======================
function S = statsTable(Tcell, labels, varName)
    n = numel(Tcell);
    N = zeros(n,1); Mean = zeros(n,1); Median = zeros(n,1); Std = zeros(n,1); P95 = zeros(n,1);
    for i=1:n
        x = double(Tcell{i}.(varName)); x = x(~isnan(x));
        N(i)=numel(x); Mean(i)=mean(x); Median(i)=median(x); Std(i)=std(x); P95(i)=prctile(x,95);
    end
    S = table(labels(:), N, Mean, Median, Std, P95, 'VariableNames', ...
        ["Trajectory","N","Mean","Median","Std","P95"]);
end

function [edges, centers] = commonBins(Tcell, binBy, numBins)
    allv = [];
    for i=1:numel(Tcell)
        allv = [allv; double(Tcell{i}.(binBy))]; %#ok<AGROW>
    end
    edges = linspace(min(allv), max(allv)+eps, numBins+1);
    centers = 0.5*(edges(1:end-1)+edges(2:end));
end

function [Mmean, Mstd] = binnedMeanStd(Tcell, binBy, edges, metricVar)
    numBins = numel(edges)-1; nRuns = numel(Tcell);
    Mmean = nan(numBins,nRuns); Mstd = nan(numBins,nRuns);
    for i=1:nRuns
        Ti = Tcell{i};
        b = discretize(double(Ti.(binBy)), edges);
        x = double(Ti.(metricVar));
        ok = ~isnan(b) & ~isnan(x);
        b=b(ok); x=x(ok);
        Mmean(:,i)=accumarray(b,x,[numBins 1],@mean,NaN);
        Mstd(:,i)=accumarray(b,x,[numBins 1],@std, NaN);
    end
end

function plotCDF(Tcell, labels, varName, xlab, ttl, savePath, saveFigs)
    figure('Color','w'); hold on; grid on;
    for i=1:numel(Tcell)
        x = sort(double(Tcell{i}.(varName)));
        x = x(~isnan(x));
        p = linspace(0,100,numel(x));
        plot(x,p,'LineWidth',1.5);
    end
    xlabel(xlab); ylabel("Percentile"); title(ttl);
    legend(labels,'Location','southeast'); hold off;
    if saveFigs, saveas(gcf, savePath); end
end
