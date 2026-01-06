clc; clear; close all;

PLOT_FIGURES = true;
NUM_REPEATS = 1;

opts.MaxIterations = 100;%5000;
opts.MinCostThreshold = 1;%1e-6;     % 1e-3
opts.RelCostDelta = 1;%1e-6;         % 1e-6
opts.StabilityPatience = 10;      % 500
opts.MaxTrial = 250;

% folderNameList = {'DataFolder_1920x1080_50_NOISY_REL_MOTION_TRANS_XYZ_PAN_3D_EIGHT2D', ...
%                   'DataFolder_1920x10850_50_NOISY_REL_MOTION_TRANS_XYZ_PAN_3D_APP_EIGHT2D', ...
%                   'DataFolder_1920x1080_50_NOISY_REL_MOTION_TRANS_XYZ_PAN_3D_EIGHT3D'};


folderName = 'D:\GitHub\PoseFromTrajectory_GeometricSolution_Matlab';
transErrList = {[0.0, 0.0, 0.0]};
angErrList = 0.0;

logFileName = PoseEstimationICIP_RealData(folderName, transErrList, angErrList, PLOT_FIGURES, NUM_REPEATS, opts);


%stats  = plotResults_updated_thetaOnly_perfRuntime(['Logs/' logFileName], false);
%figure,

%plotResults_userCorrected(['Logs/' logFileName], false);


