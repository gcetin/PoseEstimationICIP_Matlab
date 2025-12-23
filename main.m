clc; clear; close all;

PLOT_FIGURES = false;
NUM_REPEATS = 5;

folderNameList = {'DataFolder_1920x1080_100_NOISY_REL_MOTION_TRANS_XYZ_PAN_3D_EIGHT'};
transErrList = {[0.0, 0.0, 0.0]};
angErrList = 0.0;


logFileName = PoseEstimationICIP(folderNameList, transErrList, angErrList, PLOT_FIGURES, NUM_REPEATS);

plotResults(['Logs/' logFileName])


