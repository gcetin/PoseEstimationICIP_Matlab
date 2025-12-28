function logFileName = PoseEstimationICIP(folderNameList, transErrList, angErrList, PLOT_FIGURES, NUM_REPEATS, opts)
    for f = 1:length(folderNameList)
        for t = 1:length(transErrList)
            for a = 1:length(angErrList)
                logFileName = run_experiment(folderNameList{f}, angErrList(a), transErrList{t}, PLOT_FIGURES, NUM_REPEATS, opts);
            end
        end
    end
end


function logFileName = run_experiment(expFolderName, angStd, transStd, PLOT_FIGURES, NUM_REPEATS, opts)
  
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
    [gtDataDict, trajData, idealTgtTraj3d, noisyTgtTraj3d] = DataLoader.LoadInOutData(expDirName);

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
    % MAX_ITERATIONS = 5000;
    % MIN_COST_THRESHOLD = 20.0;%1e-3;
    % REL_COST_DELTA = 1e-3;%1e-6;
    % COST_STABILITY_PATIENCE = 10; %500

    % Optimization Settings
    % opts.MaxIterations = 5000;
    % opts.MinCostThreshold = 20.0;     % 1e-3
    % opts.RelCostDelta = 1e-3;         % 1e-6
    % opts.StabilityPatience = 10;      % 500
    % opts.MaxTrial = 1000;
    
    %MAX_TRIAL = 1000;

    % Logging setup (Skipped file creation for brevity, printing to console)
    fprintf('Starting Experiment: %s\n', expFolderName);

    % --- LOGGING SETUP START ---
    logFileName = UtilityFunctions.build_log_file_name(...
        expFolderName, angStd, transStd, opts.MaxIterations, opts.MinCostThreshold);
    
    fprintf('Logging to: Logs/%s\n', logFileName);
    logFid = UtilityFunctions.createLogDataFiles(logFileName);
    % --- LOGGING SETUP END ---    

    % NUM_REPEATS = 1;
    resNormArr = [];
    noisyTrajRmseList = [];


    % Levenberg-Marquardt Options
    options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', ...
                       'Display', 'off', 'MaxFunctionEvaluations', 10000, ...
                       'FunctionTolerance', 1e-12);

    for r = 1:NUM_REPEATS
        for i = 1:numSamples
            % Get observed 2D points (N x 2)
            point2d = squeeze(trajData(i, :, :));
            numPoints = size(point2d, 1);
            
            idealTraj = squeeze(idealTgtTraj3d(i,:,:));
            noisytTraj = squeeze(noisyTgtTraj3d(i,:,:));

            noisyMotionMetric = UtilityFunctions.pointwiseErrors(idealTraj, noisytTraj, false);
            noisyTrajRmse = noisyMotionMetric.rmse;

            % Precompute rays (Identity camera)
            [rayOrigins, rayDirections] = UtilityFunctions.precompute_camera_rays_identity(point2d, Params.KK);

            % Ground Truth Matrix
            gtTransMat = UtilityFunctions.getGtTransMat(gtDataDict, i);

            bestCost = Inf;
            bestCostList = [];
            badCostList = [];
            tic;

            for trialIdx=1:opts.MaxTrial 
            
                % Add Noise for initialization
                [noisyGtTransMat, yawErr, transErr] = UtilityFunctions.addNoiseToPose(gtTransMat, angStd, transStd);
                currentTransMat = noisyGtTransMat;
                %currentTransMat = gtTransMat;
    
                % % Project
                % % 1. Transform 3D points
                % pts3d_h = [idealPts3d; ones(1, size(idealPts3d,2))];
                % 
                % T = gtTransMat;
                % 
                % p_cam = T * pts3d_h;
                % p_cam = p_cam(1:3, :) ./ p_cam(4, :);
                % 
                % % 2. Project to Pixel
                % uv_hom = Params.KK * p_cam;
                % uv = uv_hom(1:2, :) ./ uv_hom(3, :);
                % 
                % figure,
                % plot(point2d(:,1), point2d(:,2), 'ob'), hold on,
                % plot(uv(1,:), uv(2,:), '*r'),
                % axis([0, 1920, 0, 1024])
                % aaaa
    
    
                % Init Loop variables
                current_cost = inf;
                prev_cost = inf;
                stability_counter = 0;
                
                % Initial Guess Vector (Rodrigues 1x3, Trans 1x3)
                initPoseGuess = UtilityFunctions.formInitialPoseGuess(currentTransMat);               
            
                % --- Iterative Refinement ---
                iteration = 0;
                while (iteration < opts.MaxIterations) && (current_cost > opts.MinCostThreshold) && (stability_counter < opts.StabilityPatience)
                    
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
                    % x_opt
                    % initPoseGuess
                    % resnorm
                    % disp("--------------------------")
                    
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
                        if rel_imp < opts.RelCostDelta
                            stability_counter = stability_counter + 1;
                        else
                            stability_counter = 0;
                        end
                    end
                    prev_cost = current_cost;
                    iteration = iteration + 1;
                end
                
                % --- CHECK IF THIS TRIAL IS THE BEST ---
                if current_cost < bestCost
                    bestCost = current_cost;
                    bestTransMat = currentTransMat;                    
                    bestCostList = [bestCostList bestCost];
                    bestNoisyTransMat = noisyGtTransMat;
                else
                    badCostList = [badCostList current_cost];
                end                
            end
            duration = toc;
            % figure(13), plot(bestCostList)
            % figure(14), plot(badCostList)
            
            %fprintf('Sample %d | Iter: %d | Cost: %.4f | StabilityCounter: %d | yawErr: %.1f | transErrX: %.1f | transErrY: %.1f | transErrZ: %.1f\n', i, iteration, bestCost, stability_counter, yawErr, transErr(1), transErr(2), transErr(3));


            %if stability_counter == 100
                %fprintf('Sample %d | Iter: %d | Cost: %.4f | StabilityCounter: %d | yawErr: %.1f | transErrX: %.1f | transErrY: %.1f | transErrZ: %.1f\n', i, iteration, bestCost, stability_counter, yawErr, transErr(1), transErr(2), transErr(3));
                %if bestCost > 100 
                    fprintf('Sample %d | Iter: %d | Cost: %.4f | StabilityCounter: %d | yawErr: %.1f | transErrX: %.1f | transErrY: %.1f | transErrZ: %.1f | noisyTrajRmse: %.2f\n', i, iteration, bestCost, stability_counter, yawErr, transErr(1), transErr(2), transErr(3), noisyTrajRmse);
                    
                    % 2. Extract Rotation Matrix (3x3) and Translation Vector (3x1)
                    R_gt = gtTransMat(1:3, 1:3);
                    t_gt = gtTransMat(1:3, 4);

                    R_noisy = bestNoisyTransMat(1:3, 1:3);
                    t_noisy = bestNoisyTransMat(1:3, 4);

                    R_est = bestTransMat(1:3, 1:3);
                    t_est = bestTransMat(1:3, 4);                    

                    % 3. Convert Rotation to Euler Angles
                    % Using 'ZYX' convention (Yaw, Pitch, Roll) to match your previous code
                    eulRad = rotm2eul(R_gt, 'ZYX'); % Returns [Yaw, Pitch, Roll] in radians
                    eulDeg = rad2deg(eulRad);       % Convert to degrees

                    eulRadN = rotm2eul(R_noisy, 'ZYX'); % Returns [Yaw, Pitch, Roll] in radians
                    eulDegN = rad2deg(eulRadN);       % Convert to degrees

                    eulRadEst = rotm2eul(R_est, 'ZYX'); % Returns [Yaw, Pitch, Roll] in radians
                    eulDegEst = rad2deg(eulRadEst);       % Convert to degrees                    

                    % Map to specific names for clarity
                    yaw_deg   = eulDeg(1); % Z-axis
                    pitch_deg = eulDeg(2); % Y-axis
                    roll_deg  = eulDeg(3); % X-axis

                    yawN_deg   = eulDegN(1); % Z-axis
                    pitchN_deg = eulDegN(2); % Y-axis
                    rollN_deg  = eulDegN(3); % X-axis

                    yawEst_deg   = eulDegEst(1); % Z-axis
                    pitchEst_deg = eulDegEst(2); % Y-axis
                    rollEst_deg  = eulDegEst(3); % X-axis                    
                    
                    % 4. Display Results
                    fprintf('--- Translation (gt - noisy - estimated - abs(error)) ---\n');
                    fprintf('X: (%.4f) - (%.4f) - (%.4f) - (%.4f)\n', t_gt(1), t_noisy(1), t_est(1), abs(t_gt(1)-t_est(1)));
                    fprintf('Y: (%.4f) - (%.4f) - (%.4f) - (%.4f)\n', t_gt(2), t_noisy(2), t_est(2), abs(t_gt(2)-t_est(2)));
                    fprintf('Z: (%.4f) - (%.4f) - (%.4f) - (%.4f)\n', t_gt(3), t_noisy(3), t_est(3), abs(t_gt(3)-t_est(3)));
                    
                    fprintf('\n--- Rotation (Degrees)  (gt - noisy - estimated - abs(error))  ---\n');
                    % fprintf('Roll  (X): (%.4f) - (%.4f) - (%.4f)\n', roll_deg, rollN_deg, rollEst_deg);
                    fprintf('Yaw (Y): (%.4f) - (%.4f) - (%.4f) - (%.4f)\n', pitch_deg, pitchN_deg, pitchEst_deg, abs(pitch_deg - pitchEst_deg));
                    % fprintf('Pitch   (Z): (%.4f) - (%.4f) - (%.4f)\n', yaw_deg, yawN_deg, yawEst_deg);    

                    fprintf('---------------------------------------------------------------------------\n');
                %end
                %end
                resNormArr = [resNormArr bestCost];
                noisyTrajRmseList = [noisyTrajRmseList noisyTrajRmse];              



            % end
            % --- FINAL VISUALIZATION ---
            if PLOT_FIGURES && bestCost > 100.0%stability_counter == opts.StabilityPatience
                % 1. Transform cloud to FINAL estimated pose
                moreDense3dInCamFrame = UtilityFunctions.convertWorldPtIntoCamFrame(moreDense3d, currentTransMat);

                % 2. Re-calculate all matches for the final pose
                % (We need the indices to draw lines to the correct points)
                all_best_indices = zeros(numPoints, 1);

                % Fast loop to get indices (calculation only, no plotting)
                for k = 1:numPoints
                    [~, idx] = UtilityFunctions.project_cloud_onto_ray_fast(rayOrigins(k,:), rayDirections(k,:), moreDense3dInCamFrame');
                    all_best_indices(k) = idx;
                end              

                 % yawErr, transErr


                % 3. Visualize Everything At Once
                fprintf('Plotting full scene for Sample %d...\n', i);
                UtilityFunctions.visualize_entire_scene(rayOrigins, rayDirections, moreDense3dInCamFrame, all_best_indices, gtTransMat, currentTransMat);
                % Optional: Pause so you can inspect it before the next sample starts
                pause;
            end


            % % --- VISUALIZATION AFTER OPTIMIZATION ---
            % if PLOT_FIGURES
            %     % 1. Update cloud to FINAL pose
            %     moreDense3dInCamFrame = UtilityFunctions.convertWorldPtIntoCamFrame(moreDense3d, currentTransMat);
            % 
            %     % 2. Initialize Figure ONCE before the loop
            %     figure(100); 
            %     clf;       % Clear any old data
            %     hold on;   % Start holding plots
            %     title(sprintf('Growing Scene: Sample %d', i));
            % 
            %     fprintf('Building scene for Sample %d...\n', i);
            % 
            %     % 3. Loop to grow the scene
            %     for k = 1:numPoints
            %         % Find best match for this ray in the final cloud
            %         [~, best_idx] = UtilityFunctions.project_cloud_onto_ray_fast(rayOrigins(k,:), rayDirections(k,:), moreDense3dInCamFrame');
            % 
            %         % Add this ray and its points to the existing plot
            %         UtilityFunctions.visualize_ray_projection(rayOrigins(k,:), rayDirections(k,:), moreDense3dInCamFrame, best_idx);
            % 
            %         % Optional: Adjust view to follow the action?
            %         % view(3); 
            % 
            %         % Fast pause for animation effect
            %         pause(0.01); 
            %     end
            %     hold off; % Stop holding after loop finishes
            % end               

            % <--- INSERT THIS LINE HERE --->
            % Visualize the projection for the last processed ray (j) and its best match
            % Note: moreDense3dInCamFrame is 3xN, matching the new function signature.
            %UtilityFunctions.visualize_ray_projection(rayOrigins(j,:), rayDirections(j,:), moreDense3dInCamFrame, best_idx);            

            % --- LOGGING DATA START ---
            
            % 1. Determine Validity flags
            converged_on_delta = (stability_counter >= opts.StabilityPatience);
            lossCond = (current_cost <= opts.MinCostThreshold) || converged_on_delta;
            ctrCond = (iteration >= opts.MaxIterations);
            
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
            R_est_final = bestTransMat(1:3, 1:3);
            t_est_final = bestTransMat(1:3, 4);

            % 4. Write to file
            UtilityFunctions.logData(logFid, (r-1)*numSamples + i, resultValidity, ...
                lmRes, gtTransMat, bestNoisyTransMat, ...
                R_est_final, t_est_final, duration, ...
                lossCond, ctrCond);
            
            % --- LOGGING DATA END ---

            
            % Logging (Brief)
            % if mod(i, 100) == 0
            %     fprintf('Sample %d | Iter: %d | Cost: %.4f | StabilityCounter: %d\n', i, iteration, current_cost, stability_counter);
            % end
        end
    end
    disp('DONE.');
    fclose(logFid);

    figure(123214), plot(resNormArr)
    figure(121234), plot(noisyTrajRmseList)

end