classdef DataLoader
    methods (Static)
        function [gtDataDict, trajData] = LoadInOutData(dirName)
            % 1. Load Ground Truth
            gtFile = fullfile(dirName, 'ground_truth_data.txt');
            opts = detectImportOptions(gtFile, 'FileType', 'text');
            gtTable = readtable(gtFile, opts);
            
            % Convert to struct
            gtDataDict.X = gtTable.ideal_distance_x_m_;
            gtDataDict.Y = gtTable.ideal_distance_y_m_;
            gtDataDict.Z = gtTable.ideal_distance_z_m_;
            gtDataDict.qw = gtTable.ideal_quat_w;
            gtDataDict.qx = gtTable.ideal_quat_x;
            gtDataDict.qy = gtTable.ideal_quat_y;
            gtDataDict.qz = gtTable.ideal_quat_z;
            
            numSamples = height(gtTable);
            
            % 2. Load Trajectories
            % Read first file to get dimensions
            temp = readmatrix(fullfile(dirName, 'TrajectoryData_1.txt'), 'FileType', 'text');
            [T, ~] = size(temp);
            
            trajData = zeros(numSamples, T, 2);
            
            % Can use parfor here for speed if many files
            parfor i = 1:numSamples
                fname = fullfile(dirName, sprintf('TrajectoryData_%d.txt', i));
                data = readmatrix(fname, 'FileType', 'text');
                % Columns 1,2 are x,y (0-based index in description, 1,2 here)
                trajData(i, :, :) = data(:, 1:2);
            end
        end
    end
end