classdef DataLoader_RealData
    methods (Static)
        function [gtDataDict, trajData, idealTgtTraj3d, noisyTgtTraj3d] = LoadInOutData(dirName)
            % 1. Load Ground Truth
            % Convert to struct
            gtDataDict.X = [0.0, 0.0, 0.0];
            gtDataDict.Y = [0.0, 0.0, 0.0];
            gtDataDict.Z = [15.0, 20.0, 25.0];
            gtDataDict.r0 = [0.0, 0.0, 0.0];
            gtDataDict.r1 = [0.0, 0.0, 0.0];
            gtDataDict.r2 = [0.0, 0.0, 0.0];    
            
            numSamples = 3;
            fileNames = {'FlightTrajData_15m.txt', 'FlightTrajData_20m.txt', 'FlightTrajData_25m.txt'};


            % 2. Load Trajectories
            % Read first file to get dimensions
            temp = readmatrix(fullfile(dirName, 'FlightTrajData_15m.txt'), 'FileType', 'text');
            [T, ~] = size(temp);
            
            trajData = zeros(numSamples, T, 2);
            
            % Can use parfor here for speed if many files
            for i = 1:numSamples
                fname = fullfile(dirName, fileNames{i});
                data = readmatrix(fname, 'FileType', 'text');
                % Columns 1,2 are x,y (0-based index in description, 1,2 here)
                trajData(i, :, :) = data(:, 1:2);
            end
        end
    end
end