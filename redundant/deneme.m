clear; clc; close all;

%% FILE
inFile = ".\Logs\DataFolder_1920x1080_250_NOISY_REL_MOTION_TRANS_XYZ_PAN_3D_EIGHT3D_5000it_1e-06minCost_250maxTrial_ang5_trans2_2_5.txt";
outDir = "analysis_outputs_matlab_single";
if ~exist(outDir,"dir"); mkdir(outDir); end

%% LOAD
opts = detectImportOptions(inFile,"Delimiter","\t");
opts = setvartype(opts,"double");
T = readtable(inFile,opts);

fprintf("Loaded %d samples\n",height(T));

%% FILTER SUCCESS
if any(strcmp(T.Properties.VariableNames,"Success"))
    T = T(T.Success==1,:);
    fprintf("Using %d successful samples\n",height(T));
end

%% ERRORS
T.dRoll  = wrapTo180(T.Est_Roll  - T.GT_Roll);
T.dPitch = wrapTo180(T.Est_Pitch - T.GT_Pitch);
T.dYaw   = wrapTo180(T.Est_Yaw   - T.GT_Yaw);

T.tErr_x = T.Est_x - T.GT_x;
T.tErr_y = T.Est_y - T.GT_y;
T.tErr_z = T.Est_z - T.GT_z;
T.tErrMag = sqrt(T.tErr_x.^2 + T.tErr_y.^2 + T.tErr_z.^2);

if ~any(strcmp(T.Properties.VariableNames,"RotErrMag"))
    T.RotErrMag = sqrt(T.dRoll.^2 + T.dPitch.^2 + T.dYaw.^2);
    warning("RotErrMag missing → Euler L2 proxy used");
end

%% AGGREGATE STATS (INLINE)
names = ["dRoll(deg)","dPitch(deg)","dYaw(deg)","RotErrMag","tErrMag"];
data  = {T.dRoll, T.dPitch, T.dYaw, T.RotErrMag, T.tErrMag};

Metric = strings(numel(names),1);
N = zeros(numel(names),1);
RMSE = zeros(numel(names),1);
Median = zeros(numel(names),1);
P90 = zeros(numel(names),1);
P95 = zeros(numel(names),1);
P99 = zeros(numel(names),1);

for i=1:numel(names)
    x = double(data{i}(:));
    Metric(i) = names(i);
    N(i) = numel(x);
    RMSE(i) = sqrt(mean(x.^2));
    Median(i) = median(x);
    P90(i) = prctile(x,90);
    P95(i) = prctile(x,95);
    P99(i) = prctile(x,99);
end

stats = table(Metric,N,RMSE,Median,P90,P95,P99);
disp(stats);
writetable(stats,fullfile(outDir,"aggregate_stats.csv"));

%% SCATTER: Noise vs Error
if all(ismember(["NoiseAngMag","RotErrMag"],T.Properties.VariableNames))
    figure; grid on;
    scatter(T.NoiseAngMag,T.RotErrMag,10,'filled','MarkerFaceAlpha',0.35);
    xlabel("NoiseAngMag (deg)");
    ylabel("RotErrMag");
    title("Noise angle vs rotation error");
    saveas(gcf,fullfile(outDir,"noiseAng_vs_rotErr.png"));
end

if all(ismember(["NoiseTransMag","tErrMag"],T.Properties.VariableNames))
    figure; grid on;
    scatter(T.NoiseTransMag,T.tErrMag,10,'filled','MarkerFaceAlpha',0.35);
    xlabel("NoiseTransMag");
    ylabel("Translation error magnitude");
    title("Noise translation vs translation error");
    saveas(gcf,fullfile(outDir,"noiseTrans_vs_tErr.png"));
end

%% PERCENTILE (CDF) PLOTS
x = sort(T.RotErrMag);
p = linspace(0,100,numel(x));
figure; grid on;
plot(x,p,'LineWidth',1.5);
xlabel("RotErrMag");
ylabel("Percentile");
title("Rotation error distribution");
saveas(gcf,fullfile(outDir,"rotErr_percentile.png"));

x = sort(T.tErrMag);
p = linspace(0,100,numel(x));
figure; grid on;
plot(x,p,'LineWidth',1.5);
xlabel("tErrMag");
ylabel("Percentile");
title("Translation error distribution");
saveas(gcf,fullfile(outDir,"tErr_percentile.png"));

%% HISTOGRAMS
figure; grid on;
histogram(T.RotErrMag,60);
xlabel("RotErrMag"); ylabel("Count");
title("Histogram: rotation error");
saveas(gcf,fullfile(outDir,"hist_rotErr.png"));

figure; grid on;
histogram(T.tErrMag,60);
xlabel("tErrMag"); ylabel("Count");
title("Histogram: translation error");
saveas(gcf,fullfile(outDir,"hist_tErr.png"));

fprintf("✔ Analysis finished. Results saved to %s\n",outDir);

%% =========================================================
%% GT VALUES vs ERROR VALUES (BINNED BAR PLOTS) - FIXED
%% =========================================================

numBins = 20;

figure('Color','w','Name','GT Translation vs Error (Binned)');

% ---------- X ----------
subplot(3,1,1);
edges = linspace(min(T.GT_x), max(T.GT_x)+eps, numBins+1);
[~,~,binIdx] = histcounts(T.GT_x, edges);

errMean = accumarray(binIdx(binIdx>0), abs(T.tErr_x(binIdx>0)), ...
                     [numBins 1], @mean, NaN);

bar(errMean, 'FaceAlpha',0.7);   % <<< NO X VALUES
grid on;
xlabel('GT X bins');
ylabel('Mean |Error X|');
title('GT X vs Translation Error');
xticks(1:numBins);
xticklabels(round(linspace(min(T.GT_x), max(T.GT_x), numBins),2));
xtickangle(45);

% ---------- Y ----------
subplot(3,1,2);
edges = linspace(min(T.GT_y), max(T.GT_y)+eps, numBins+1);
[~,~,binIdx] = histcounts(T.GT_y, edges);

errMean = accumarray(binIdx(binIdx>0), abs(T.tErr_y(binIdx>0)), ...
                     [numBins 1], @mean, NaN);

bar(errMean, 'FaceAlpha',0.7);
grid on;
xlabel('GT Y bins');
ylabel('Mean |Error Y|');
title('GT Y vs Translation Error');
xticks(1:numBins);
xticklabels(round(linspace(min(T.GT_y), max(T.GT_y), numBins),2));
xtickangle(45);

% ---------- Z ----------
subplot(3,1,3);
edges = linspace(min(T.GT_z), max(T.GT_z)+eps, numBins+1);
[~,~,binIdx] = histcounts(T.GT_z, edges);

errMean = accumarray(binIdx(binIdx>0), abs(T.tErr_z(binIdx>0)), ...
                     [numBins 1], @mean, NaN);

bar(errMean, 'FaceAlpha',0.7);
grid on;
xlabel('GT Z bins');
ylabel('Mean |Error Z|');
title('GT Z vs Translation Error');
xticks(1:numBins);
xticklabels(round(linspace(min(T.GT_z), max(T.GT_z), numBins),2));
xtickangle(45);

saveas(gcf, fullfile(outDir,'GT_vs_Translation_Error_binned_bar.png'));
