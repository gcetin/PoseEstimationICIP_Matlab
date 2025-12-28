classdef DataLoader
    methods (Static)
        function [gtDataDict, trajData, idealTgtTraj3d, noisyTgtTraj3d] = LoadInOutData(dirName)
            % 1. Load Ground Truth
            gtFile = fullfile(dirName, 'ground_truth_data.txt');
            opts = detectImportOptions(gtFile, 'FileType', 'text');
            gtTable = readtable(gtFile, opts);
            
            % Convert to struct
            gtDataDict.X = gtTable.ideal_distance_x_m_;
            gtDataDict.Y = gtTable.ideal_distance_y_m_;
            gtDataDict.Z = gtTable.ideal_distance_z_m_;
            gtDataDict.r0 = gtTable.ideal_rodrigues_x_rad_;
            gtDataDict.r1 = gtTable.ideal_rodrigues_y_rad_;
            gtDataDict.r2 = gtTable.ideal_rodrigues_z_rad_;    
            
            numSamples = height(gtTable);
            
            % 2. Load Trajectories
            % Read first file to get dimensions
            temp = readmatrix(fullfile(dirName, 'TrajectoryData_1.txt'), 'FileType', 'text');
            [T, ~] = size(temp);
            
            trajData = zeros(numSamples, T, 2);
            idealTgtTraj3d = zeros(numSamples, T, 3);
            noisyTgtTraj3d = zeros(numSamples, T, 3);
            
            % Can use parfor here for speed if many files
            parfor i = 1:numSamples
                fname = fullfile(dirName, sprintf('TrajectoryData_%d.txt', i));
                data = readmatrix(fname, 'FileType', 'text');
                % Columns 1,2 are x,y (0-based index in description, 1,2 here)
                trajData(i, :, :) = data(:, 1:2);
                idealTgtTraj3d(i, :, :) = data(:, 6:8);
                noisyTgtTraj3d(i, :, :) = data(:, 9:11);
            end
        end
    end
end