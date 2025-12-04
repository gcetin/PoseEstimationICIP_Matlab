function PoseEstimationICIP()
    % Main entry point for Pose Estimation
    clc; clear; close all;

    % Configuration
    folderNameList = {'DataFolder_1920x1080_500_IDEAL_REL_MOTION_TRANS_XYZ_PAN_3D_EIGHT'};
    transErrList = {[1.0, 1.0, 5.0]};
    angErrList = [0.0];

    for f = 1:length(folderNameList)
        for t = 1:length(transErrList)
            for a = 1:length(angErrList)
                run_experiment(folderNameList{f}, angErrList(a), transErrList{t});
            end
        end
    end
end


function run_experiment(expFolderName, angStd, transStd)
    PLOT_FIGURES = true;
    % Setup Paths
    PARTITION_NAME = 'D:\';
    BASE_DIR_NAME = 'DOKTORA\EXPERIMENTS\Dataset_Text_Python\';
    expDirName = fullfile(PARTITION_NAME, BASE_DIR_NAME, expFolderName, 'Train');

    if ~exist(expDirName, 'dir')
        error('Given folder does not exist: %s', expDirName);
    end

    % Constants
    Params = GenSyntheticDataTypes();
    
    % Load Data
    [gtDataDict, trajData] = DataLoader.LoadInOutData(expDirName);
    
    % Parse config from folder name
    [numSamples, ~, appendFlag, noiseFlag, shapeType] = UtilityFunctions.parse_folder_name(expFolderName);
    
    % Trajectory Generation (Ideal)
    deltaT = Params.deltaT;
    if appendFlag; deltaT = deltaT * 2; end
    
    % Create dense trajectory
    moreDense3d = TrajectoryGenerator.create_dense_trajectory(shapeType, appendFlag, deltaT);
    
    % Sample ideal points (downsample ratio 5)
    idealPts3d = TrajectoryGenerator.sample_trajectory(moreDense3d, 5);

    % Optimization Settings
    MAX_ITERATIONS = 200;
    MIN_COST_THRESHOLD = 1000.0;
    REL_COST_DELTA = 1e-3;
    COST_STABILITY_PATIENCE = 5;

    % Logging setup (Skipped file creation for brevity, printing to console)
    fprintf('Starting Experiment: %s\n', expFolderName);

    % --- LOGGING SETUP START ---
    logFileName = UtilityFunctions.build_log_file_name(...
        expFolderName, angStd, transStd, MAX_ITERATIONS, MIN_COST_THRESHOLD);
    
    fprintf('Logging to: Logs/%s\n', logFileName);
    logFid = UtilityFunctions.createLogDataFiles(logFileName);
    % --- LOGGING SETUP END ---    

    NUM_REPEATS = 1;

    for r = 1:NUM_REPEATS
        for i = 1:numSamples
            % Get observed 2D points (N x 2)
            point2d = squeeze(trajData(i, :, :));
            numPoints = size(point2d, 1);

            % Precompute rays (Identity camera)
            [rayOrigins, rayDirections] = UtilityFunctions.precompute_camera_rays_identity(point2d, Params.KK);

            % Ground Truth Matrix
            gtTransMat = UtilityFunctions.getGtTransMat(gtDataDict, i);
            
            % Add Noise for initialization
            noisyGtTransMat = UtilityFunctions.addNoiseToPose(gtTransMat, angStd, transStd);
            currentTransMat = noisyGtTransMat;

            % Init Loop variables
            current_cost = inf;
            prev_cost = inf;
            stability_counter = 0;
            
            % Initial Guess Vector (Rodrigues 1x3, Trans 1x3)
            initPoseGuess = UtilityFunctions.formInitialPoseGuess(currentTransMat);
            
            % Levenberg-Marquardt Options
            options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', ...
                                   'Display', 'off', 'MaxFunctionEvaluations', 100000, ...
                                   'FunctionTolerance', 1e-8);
            tic;
            % --- Iterative Refinement ---
            iteration = 0;
            while (iteration < MAX_ITERATIONS) && (current_cost > MIN_COST_THRESHOLD) && (stability_counter < COST_STABILITY_PATIENCE)
                
                % 1. Convert Dense Cloud to Camera Frame using CURRENT Pose
                % (3 x N)
                moreDense3dInCamFrame = UtilityFunctions.convertWorldPtIntoCamFrame(moreDense3d, currentTransMat);
                
                % 2. Correspondence Establishment
                % For every ray, find the closest point in the dense cloud
                closestPointToTheRayCoord = zeros(3, numPoints);
                
                % Vectorized Ray-Cloud Intersection
                for j = 1:numPoints
                    [~, best_idx] = UtilityFunctions.project_cloud_onto_ray_fast(rayOrigins(j,:), rayDirections(j,:), moreDense3dInCamFrame');
                    
                    closestPtCam = moreDense3dInCamFrame(:, best_idx);
                    
                    % Transform back to Object Frame
                    closestPointToTheRayCoord(:, j) = UtilityFunctions.transform_cam_point_to_object(currentTransMat, closestPtCam);
                end

                % 3. Optimize Pose
                % We optimize the delta or the absolute? Python code optimizes 'initPoseGuess' 
                % against the NEW correspondences.
                
                fun = @(x) UtilityFunctions.residuals(x, closestPointToTheRayCoord, point2d, Params.KK, Params.distCoeffs);
                
                [x_opt, resnorm, ~, exitflag] = lsqnonlin(fun, initPoseGuess, [], [], options);
                
                current_cost = resnorm; % Squared error sum

                % Update Pose (Python code extracts result.x)
                optRotVec = x_opt(1:3);
                optTransVec = x_opt(4:6);
                
                % Update Matrix (Optimization returns absolute pose for the correspondences)
                R_opt = UtilityFunctions.rodrigues(optRotVec);
                currentTransMat = eye(4);
                currentTransMat(1:3, 1:3) = R_opt;
                currentTransMat(1:3, 4) = optTransVec';
                
                % Update Guess for next iteration
                initPoseGuess = x_opt;

                % Convergence Check
                if isfinite(prev_cost)
                    rel_imp = abs(prev_cost - current_cost) / max(abs(prev_cost), 1e-12);
                    if rel_imp < REL_COST_DELTA
                        stability_counter = stability_counter + 1;
                    else
                        stability_counter = 0;
                    end
                end
                prev_cost = current_cost;
                iteration = iteration + 1;
            end
            toc;


            % --- VISUALIZATION AFTER OPTIMIZATION ---
            if PLOT_FIGURES
                % 1. Update cloud to FINAL pose
                moreDense3dInCamFrame = UtilityFunctions.convertWorldPtIntoCamFrame(moreDense3d, currentTransMat);
                
                % 2. Initialize Figure ONCE before the loop
                figure(100); 
                clf;       % Clear any old data
                hold on;   % Start holding plots
                title(sprintf('Growing Scene: Sample %d', i));
                
                fprintf('Building scene for Sample %d...\n', i);
                
                % 3. Loop to grow the scene
                for k = 1:numPoints
                    % Find best match for this ray in the final cloud
                    [~, best_idx] = UtilityFunctions.project_cloud_onto_ray_fast(rayOrigins(k,:), rayDirections(k,:), moreDense3dInCamFrame');
                    
                    % Add this ray and its points to the existing plot
                    UtilityFunctions.visualize_ray_projection(rayOrigins(k,:), rayDirections(k,:), moreDense3dInCamFrame, best_idx);
                    
                    % Optional: Adjust view to follow the action?
                    % view(3); 
                    
                    % Fast pause for animation effect
                    pause(0.01); 
                end
                hold off; % Stop holding after loop finishes
            end               

            % <--- INSERT THIS LINE HERE --->
            % Visualize the projection for the last processed ray (j) and its best match
            % Note: moreDense3dInCamFrame is 3xN, matching the new function signature.
            UtilityFunctions.visualize_ray_projection(rayOrigins(j,:), rayDirections(j,:), moreDense3dInCamFrame, best_idx);            

            % --- LOGGING DATA START ---
            
            % 1. Determine Validity flags
            converged_on_delta = (stability_counter >= COST_STABILITY_PATIENCE);
            lossCond = (current_cost <= MIN_COST_THRESHOLD) || converged_on_delta;
            ctrCond = (iteration >= MAX_ITERATIONS);
            
            resultValidity = 0;
            if lossCond
                resultValidity = 1;
            end
            
            % 2. Pack LM Results
            % lsqnonlin returns resnorm (cost), output struct (iterations, funcCount)
            % 'output' variable comes from lsqnonlin: [x,resnorm,residual,exitflag,output]
            [~, resnorm, ~, ~, output] = lsqnonlin(fun, initPoseGuess, [], [], options);
            
            lmRes.cost = resnorm;
            lmRes.optimality = output.firstorderopt;
            lmRes.funcCount = output.funcCount;
            lmRes.iterations = output.iterations;

            % 3. Extract final Estimated Rotation/Translation from Matrix
            % (Assuming currentTransMat was updated with final result)
            R_est_final = currentTransMat(1:3, 1:3);
            t_est_final = currentTransMat(1:3, 4);

            % 4. Write to file
            UtilityFunctions.logData(logFid, (r-1)*numSamples + i, resultValidity, ...
                lmRes, gtTransMat, noisyGtTransMat, ...
                R_est_final, t_est_final, duration, ...
                lossCond, ctrCond);
            
            % --- LOGGING DATA END ---

            
            % Logging (Brief)
            if mod(i, 100) == 0
                fprintf('Sample %d | Iter: %d | Cost: %.4f\n', i, iteration, current_cost);
            end
        end
    end
    disp('DONE.');
end