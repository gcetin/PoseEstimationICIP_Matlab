clc; clear; close all;

PLOT_FIGURES = false;
NUM_REPEATS = 1;

opts.MaxIterations = 5000;
opts.MinCostThreshold = 1e-6;     % 1e-3
opts.RelCostDelta = 1e-6;         % 1e-6
opts.StabilityPatience = 10;      % 500
opts.MaxTrial = 250;

folderNameList = {'DataFolder_1920x1080_250_NOISY_REL_MOTION_TRANS_XYZ_PAN_3D_EIGHT3D'};
transErrList = {[2.0, 2.0, 5.0]};
angErrList = 5.0;

logFileName = PoseEstimationICIP(folderNameList, transErrList, angErrList, PLOT_FIGURES, NUM_REPEATS, opts);


%stats  = plotResults_updated_thetaOnly_perfRuntime(['Logs/' logFileName], false);
%figure,

%plotResults_userCorrected(['Logs/' logFileName], false);


