clc; clear; close all;

PLOT_FIGURES = true;

folderNameList = {'DataFolder_1920x1080_500_IDEAL_REL_MOTION_PAN_3D_EIGHT'};
transErrList = {[1.0, 1.0, 5.0]};
angErrList = [1.0];


logFileName = PoseEstimationICIP(folderNameList, transErrList, angErrList, PLOT_FIGURES);

plotResults(['Logs/' logFileName])


