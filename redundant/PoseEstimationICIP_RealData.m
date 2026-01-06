function logFileName = PoseEstimationICIP_RealData(folderName, transErrList, angErrList, PLOT_FIGURES, NUM_REPEATS, opts)
% PoseEstimationICIP_RealData
% - NUM_REPEATS is forced to 1 (real data: process each sample once)
% - For each folder and noise setting, runs run_experiment and returns last logFileName.

    %#ok<INUSD> % PLOT_FIGURES may be unused in this snippet, kept for interface consistency

    logFileName = '';

    % Force real-data protocol: one pass per sample
    NUM_REPEATS = 1;

    for f = 1:length(folderName)
        for t = 1:length(transErrList)
            for a = 1:length(angErrList)
                logFileName = run_experiment(folderName, angErrList(a), transErrList{t}, PLOT_FIGURES, NUM_REPEATS, opts);
            end
        end
    end
end


function logFileName = run_experiment(expFolderName, angStd, transStd, PLOT_FIGURES, NUM_REPEATS, opts)
% run_experiment
% Real-data run:
%   - loops samples as i = 1:numSamples
%   - multiple trials per sample (opts.MaxTrial), keep the best
%   - logs results

    %#ok<INUSD> % PLOT_FIGURES kept for interface consistency
    NUM_REPEATS = 1; % enforce

    if ~exist(expFolderName, 'dir')
        error('Given folder does not exist: %s', expFolderName);
    end

    Params = GenRealDataTypes();
    [gtDataDict, trajData] = DataLoader_RealData.LoadInOutData(expFolderName);

    % Determine number of samples from data (no hard-code)
    numSamples = size(trajData, 1);
    if numSamples < 1
        error('trajData has no samples (size(trajData,1)=0).');
    end

    appendFlag = false;
    shapeType  = 0;

    % Generate Reference Dense Cloud once
    deltaT    = Params.deltaT;
    moreDense3d = TrajectoryGenerator.create_dense_trajectory(shapeType, appendFlag, deltaT);

    logFileName = UtilityFunctions.build_log_file_name( ...
        'realDataResults', angStd, transStd, ...
        opts.MaxIterations, opts.MinCostThreshold, opts.MaxTrial);

    resultsStore = cell(numSamples, 1);

    % Optimization Options
    options = optimoptions('lsqnonlin', ...
        'Algorithm', 'levenberg-marquardt', ...
        'Display', 'off', ...
        'MaxFunctionEvaluations', 10000, ...
        'FunctionTolerance', 1e-12);

    fprintf('\nStarting RealData experiment: %s\n', expFolderName);
    fprintf('Samples: %d | MaxTrial: %d | MaxIter(outer): %d\n', numSamples, opts.MaxTrial, opts.MaxIterations);

    % --- MAIN LOOP: one pass per sample (NUM_REPEATS=1) ---
    for i = 1:numSamples

        % Extract data for this sample
        point2d = squeeze(trajData(i, :, :));     % (N x 2) or (N x ?)
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

            iteration = 0; % outer refinement iteration
            local_output = struct('firstorderopt', 0, 'funcCount', 0, 'iterations', 0);

            while (iteration < opts.MaxIterations) && ...
                  (current_cost > opts.MinCostThreshold) && ...
                  (stability_counter < opts.StabilityPatience)

                % Transform dense cloud into camera frame
                moreDense3dInCamFrame = UtilityFunctions.convertWorldPtIntoCamFrame(moreDense3d, currentTransMat);

                % ---- Ensure cloud is in Mx3 format for project_cloud_onto_ray_fast ----
                % It might come as 3xM or Mx3 depending on your implementation.
                if size(moreDense3dInCamFrame, 2) == 3
                    cloudCam = moreDense3dInCamFrame;      % Mx3
                elseif size(moreDense3dInCamFrame, 1) == 3
                    cloudCam = moreDense3dInCamFrame.';    % Mx3
                else
                    error('Unexpected moreDense3dInCamFrame size: [%d x %d]. Expected Mx3 or 3xM.', ...
                        size(moreDense3dInCamFrame,1), size(moreDense3dInCamFrame,2));
                end

                % Correspondence: for each ray, find closest point in cloud
                closestPointToTheRayCoord = zeros(3, numPoints);

                for j = 1:numPoints
                    [~, best_idx] = UtilityFunctions.project_cloud_onto_ray_fast( ...
                        rayOrigins(j,:), rayDirections(j,:), cloudCam);

                    closestPtCam = cloudCam(best_idx, :).'; % 3x1 in camera frame

                    % Convert camera-frame point back to object/world coords as your pipeline expects
                    closestPointToTheRayCoord(:, j) = UtilityFunctions.transform_cam_point_to_object(currentTransMat, closestPtCam);
                end

                % Optimize with LM (lsqnonlin)
                fun = @(x) UtilityFunctions.residuals( ...
                    x, closestPointToTheRayCoord, point2d, Params.KK, Params.distCoeffs);

                [x_opt, resnorm, ~, ~, output] = lsqnonlin(fun, initPoseGuess, [], [], options);

                current_cost = resnorm;
                local_output = output;

                % Update pose from optimized parameters
                R_opt = UtilityFunctions.rodrigues(x_opt(1:3));
                currentTransMat = eye(4);
                currentTransMat(1:3, 1:3) = R_opt;
                currentTransMat(1:3, 4)   = x_opt(4:6).';

                initPoseGuess = x_opt;

                % Stability check
                if isfinite(prev_cost)
                    relDelta = abs(prev_cost - current_cost) / max(abs(prev_cost), 1e-12);
                    if relDelta < opts.RelCostDelta
                        stability_counter = stability_counter + 1;
                    else
                        stability_counter = 0;
                    end
                end
                prev_cost = current_cost;

                iteration = iteration + 1;
            end

            % Keep best trial for this sample
            if current_cost < bestLmRes.cost
                bestLmRes.cost        = current_cost;
                bestLmRes.optimality  = local_output.firstorderopt;
                bestLmRes.funcCount   = local_output.funcCount;
                bestLmRes.iterations  = local_output.iterations; % lsqnonlin iterations

                bestTransMat      = currentTransMat;
                bestNoisyTransMat = noisyGtTransMat;
            end
        end

        duration = toc(trialStart);

        % Pack result
        s = struct();
        s.sampleIdx        = i;
        s.lmRes            = bestLmRes;
        s.bestTransMat     = bestTransMat;
        s.bestNoisyTransMat= bestNoisyTransMat;
        s.duration         = duration;
        s.gtTransMat       = gtTransMat;
        s.converged        = (bestLmRes.cost <= opts.MinCostThreshold);

        resultsStore{i} = s;

        fprintf('Sample %d/%d done | best cost: %.6g | converged: %d | time: %.2fs\n', ...
            i, numSamples, bestLmRes.cost, s.converged, duration);
    end

    % --- Backup and Text Logging ---
    if ~exist('Logs', 'dir'); mkdir('Logs'); end
    save(fullfile('Logs', [logFileName(1:end-4), '_backup.mat']), 'resultsStore', 'opts');

    logFid = UtilityFunctions.createLogDataFiles(logFileName);

    for i = 1:numSamples
        res = resultsStore{i};

        R_est = res.bestTransMat(1:3, 1:3);
        t_est = res.bestTransMat(1:3, 4);

        % NOTE: ensure logData signature matches this call
        UtilityFunctions.logData( ...
            logFid, res.sampleIdx, res.converged, ...
            res.lmRes, res.gtTransMat, res.bestNoisyTransMat, ...
            R_est, t_est, res.duration, res.converged, false);
    end

    fclose(logFid);
    fprintf('\nExperiment %s Done.\n', expFolderName);
end
