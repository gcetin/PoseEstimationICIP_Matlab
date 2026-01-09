clear all,
clc,
close all,

logFileName = 'DataFolder_1920x1080_250_NOISY_REL_MOTION_TRANS_XYZ_PAN_3D_EIGHT3D_5000it_1e-06minCost_250maxTrial_ang5_trans2_2_5.txt';

plotResults_updated_thetaOnly_perfRuntime(['Logs/' logFileName], false);

%plotResults_userCorrected(['Logs/' logFileName], false);
