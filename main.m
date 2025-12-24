clc; clear; close all;

PLOT_FIGURES = true;
NUM_REPEATS = 1;

folderNameList = {'DataFolder_1920x1080_50_NOISY_REL_MOTION_TRANS_XYZ_PAN_3D_EIGHT'};
transErrList = {[1.0, 1.0, 3.0]};
angErrList = 3.0;


logFileName = PoseEstimationICIP(folderNameList, transErrList, angErrList, PLOT_FIGURES, NUM_REPEATS);

plotResults(['Logs/' logFileName])


