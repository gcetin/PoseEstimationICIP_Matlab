function logFileName = PoseEstimationICIP(folderNameList, transErrList, angErrList, PLOT_FIGURES, NUM_REPEATS, opts)
    logFileName = ''; 
    for f = 1:length(folderNameList)
        for t = 1:length(transErrList)
            for a = 1:length(angErrList)
                logFileName = run_experiment(folderNameList{f}, angErrList(a), transErrList{t}, PLOT_FIGURES, NUM_REPEATS, opts);
            end
        end
    end
end

function logFileName = run_experiment(expFolderName, angStd, transStd, PLOT_FIGURES, NUM_REPEATS, opts)
    % --- 1. Setup and Data Loading ---
    PARTITION_NAME = 'D:\';
    BASE_DIR_NAME = 'DOKTORA\EXPERIMENTS\Dataset_Text_Python\';
    expDirName = fullfile(PARTITION_NAME, BASE_DIR_NAME, expFolderName, 'Train');
    
    if ~exist(expDirName, 'dir')
        error('Given folder does not exist: %s', expDirName);
    end

    Params = GenSyntheticDataTypes();
    [gtDataDict, trajData, idealTgtTraj3d, noisyTgtTraj3d] = DataLoader.LoadInOutData(expDirName);

    % for i=1: size(idealTgtTraj3d,1)
    %     % Extract trajectories
    %     X = [idealTgtTraj3d(i,:,1), noisyTgtTraj3d(i,:,1)];
    %     Y = [idealTgtTraj3d(i,:,2), noisyTgtTraj3d(i,:,2)];
    %     Z = [idealTgtTraj3d(i,:,3), noisyTgtTraj3d(i,:,3)];
    % 
    %     xlim_all = [min(X), max(X)];
    %     ylim_all = [min(Y), max(Y)];
    %     %zlim_all = [min(Z), max(Z)];
    %     zlim_all = [-5, 5];
    % 
    %     figure(1),
    % 
    %     subplot(1,2,1)
    %     plot3(idealTgtTraj3d(i,:,1), idealTgtTraj3d(i,:,2), idealTgtTraj3d(i,:,3), 'b')
    %     grid on
    %     axis equal
    %     xlim(xlim_all); ylim(ylim_all); zlim(zlim_all)
    %     title('Ideal Trajectory')
    % 
    %     subplot(1,2,2)
    %     plot3(noisyTgtTraj3d(i,:,1), noisyTgtTraj3d(i,:,2), noisyTgtTraj3d(i,:,3), 'r')
    %     grid on
    %     axis equal
    %     xlim(xlim_all); ylim(ylim_all); zlim(zlim_all)
    %     title('Noisy Trajectory')
    % 
    %     pause,
    % 
    % end


    [numSamples, ~, appendFlag, ~, shapeType] = UtilityFunctions.parse_folder_name(expFolderName);
    
    % Slicing variables for parfor
    trajData = trajData;
    gtDataDict = gtDataDict; 
    
    % Generate Reference Dense Cloud once
    deltaT = Params.deltaT;
    moreDense3d = TrajectoryGenerator.create_dense_trajectory(shapeType, appendFlag, deltaT);
    
    logFileName = UtilityFunctions.build_log_file_name(expFolderName, angStd, transStd, opts.MaxIterations, opts.MinCostThreshold, opts.MaxTrial);
    
    % --- 2. Parallel Progress & ETA Setup ---
    if isempty(gcp('nocreate')); parpool; end
    
    totalIterations = NUM_REPEATS * numSamples;
    resultsStore = cell(totalIterations, 1);
    
    q = parallel.pool.DataQueue();
    numCompleted = 0;
    startTime = tic; 

    afterEach(q, @(varargin) updateProgress());
    
    function updateProgress()
        numCompleted = numCompleted + 1;
        elapsedTime = toc(startTime);
        avgTimePerTask = elapsedTime / numCompleted;
        etaSeconds = avgTimePerTask * (totalIterations - numCompleted);
        
        etaStr = datestr(etaSeconds/86400, 'HH:MM:SS');
        if mod(numCompleted, 5) == 0 || numCompleted == totalIterations
            fprintf('--> [%s] Progress: %d/%d (%.1f%%) | ETA: %s\n', ...
                datestr(now, 'HH:MM:SS'), numCompleted, totalIterations, ...
                (numCompleted/totalIterations)*100, etaStr);
        end
    end

    % Optimization Options
    options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', ...
                       'Display', 'off', 'MaxFunctionEvaluations', 10000, ...
                       'FunctionTolerance', 1e-12);

    fprintf('\nStarting Parallel experiment: %s\n', expFolderName);

    % --- 3. Parallel Loop (Flattened) ---
    parfor k = 1:totalIterations
        % Map linear index k to repeat index r and sample index i
        r = ceil(k / numSamples);
        i = mod(k-1, numSamples) + 1;
        
        % Extract data for this sample
        point2d = squeeze(trajData(i, :, :));
        numPoints = size(point2d, 1);
        [rayOrigins, rayDirections] = UtilityFunctions.precompute_camera_rays_identity(point2d, Params.KK);
        gtTransMat = UtilityFunctions.getGtTransMat(gtDataDict, i);
        
        % Trial trackers
        bestLmRes = struct('cost', Inf, 'optimality', 0, 'funcCount', 0, 'iterations', 0);
        bestTransMat = eye(4);
        bestNoisyTransMat = eye(4);
        
        trialStart = tic;
        for trialIdx = 1:opts.MaxTrial 
            [noisyGtTransMat, ~, ~] = UtilityFunctions.addNoiseToPose(gtTransMat, angStd, transStd);
            currentTransMat = noisyGtTransMat;
            current_cost = inf;
            prev_cost = inf;
            stability_counter = 0;
            initPoseGuess = UtilityFunctions.formInitialPoseGuess(currentTransMat);               
            
            iteration = 0;
            local_output = struct('firstorderopt', 0, 'funcCount', 0, 'iterations', 0);

            while (iteration < opts.MaxIterations) && (current_cost > opts.MinCostThreshold) && (stability_counter < opts.StabilityPatience)
                
                % Transform cloud
                moreDense3dInCamFrame = UtilityFunctions.convertWorldPtIntoCamFrame(moreDense3d, currentTransMat);
                
                % Correspondence
                closestPointToTheRayCoord = zeros(3, numPoints);
                for j = 1:numPoints
                    [~, best_idx] = UtilityFunctions.project_cloud_onto_ray_fast(rayOrigins(j,:), rayDirections(j,:), moreDense3dInCamFrame');
                    closestPtCam = moreDense3dInCamFrame(:, best_idx);
                    closestPointToTheRayCoord(:, j) = UtilityFunctions.transform_cam_point_to_object(currentTransMat, closestPtCam);
                end

                % Optimize
                fun = @(x) UtilityFunctions.residuals(x, closestPointToTheRayCoord, point2d, Params.KK, Params.distCoeffs);
                [x_opt, resnorm, ~, ~, output] = lsqnonlin(fun, initPoseGuess, [], [], options);
               
                current_cost = resnorm;
                local_output = output;

                % Update state
                R_opt = UtilityFunctions.rodrigues(x_opt(1:3));
                currentTransMat = eye(4);
                currentTransMat(1:3, 1:3) = R_opt;
                currentTransMat(1:3, 4) = x_opt(4:6)';
                initPoseGuess = x_opt;

                if isfinite(prev_cost)
                    if (abs(prev_cost - current_cost) / max(abs(prev_cost), 1e-12)) < opts.RelCostDelta
                        stability_counter = stability_counter + 1;
                    else
                        stability_counter = 0;
                    end
                end
                prev_cost = current_cost;
                iteration = iteration + 1;
            end
            
            % Check if this is the best trial for this sample
            if current_cost < bestLmRes.cost
                bestLmRes.cost = current_cost;
                bestLmRes.optimality = local_output.firstorderopt;
                bestLmRes.funcCount = local_output.funcCount;
                bestLmRes.iterations = local_output.iterations;
                
                bestTransMat = currentTransMat;                    
                bestNoisyTransMat = noisyGtTransMat;
            end                
        end
        duration = toc(trialStart);

        % Pack result
        s = struct();
        s.sampleIdx = i;
        s.lmRes = bestLmRes; 
        s.bestTransMat = bestTransMat;
        s.bestNoisyTransMat = bestNoisyTransMat;
        s.duration = duration;
        s.gtTransMat = gtTransMat;
        s.converged = (bestLmRes.cost <= opts.MinCostThreshold);
        
        resultsStore{k} = s;
        send(q, k);
    end

    % --- 4. Backup and Text Logging ---
    if ~exist('Logs', 'dir'); mkdir('Logs'); end
    save(fullfile('Logs', [logFileName(1:end-4), '_backup.mat']), 'resultsStore', 'opts');

    logFid = UtilityFunctions.createLogDataFiles(logFileName);
    for k = 1:totalIterations
        res = resultsStore{k};
        R_est = res.bestTransMat(1:3, 1:3);
        t_est = res.bestTransMat(1:3, 4);
        
        UtilityFunctions.logData(logFid, res.sampleIdx, res.converged, ...
            res.lmRes, res.gtTransMat, res.bestNoisyTransMat, ...
            R_est, t_est, res.duration, res.converged, false);
    end
    fclose(logFid);
    fprintf('\nExperiment %s Done.\n', expFolderName);
end