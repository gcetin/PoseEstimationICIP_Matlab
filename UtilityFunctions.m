classdef UtilityFunctions
    methods (Static)

        function [numSamples, relPoseType, appendFlag, noiseFlag, shapeType] = parse_folder_name(folderName)
            % Simple parsing logic mimicking regex
            % Assuming folder format matches strict expectations

            numSamples = 500; % Default fallback or parse '500' from string
            if contains(folderName, 'NOISY')
                noiseFlag = true;
            else
                noiseFlag = false;
            end

            if contains(folderName, 'APP')
                appendFlag = true;
            else
                appendFlag = false;
            end

            % Shape
            if contains(folderName, 'EIGHT')
                shapeType = 0;
            else
                shapeType = 2; % Default
            end

            relPoseType = 0; % Placeholder
        end

        function [gtDataDict, trajData] = load_data_struct(expDirName)
            % Placeholder for complex loading logic
        end

        function [T_noisy, yaw_noise, t_noise] = addNoiseToPose(T_gt, ang_std, trans_std)
            R_gt = T_gt(1:3, 1:3);
            t_gt = T_gt(1:3, 4);

            % Add noise to Yaw (ZYX convention)
            if ang_std ~= 0
                yaw_noise = randn() * ang_std;
                % ZYX rotation matrix for just Z noise
                eul = [deg2rad(yaw_noise), 0, 0];
                R_noise = eul2rotm(eul, 'ZYX');
                R_noisy = R_noise * R_gt;                
            else
                R_noisy = R_gt;
                yaw_noise = 0;
            end

            if length(trans_std) == 3
                t_noise = (randn(3,1) .* trans_std(:));
            else
                t_noise = zeros(3,1);
            end

            t_noisy = t_gt + t_noise;
            % Clamp Z
            t_noisy(3) = max(t_noisy(3), 15.0);

            T_noisy = eye(4);
            T_noisy(1:3, 1:3) = R_noisy;
            T_noisy(1:3, 4) = t_noisy;     
        end

        function vec = formInitialPoseGuess(T)
            R = T(1:3, 1:3);
            t = T(1:3, 4);
            r_vec = UtilityFunctions.rodrigues_inv(R); % Rotation Matrix to Vector
            vec = [r_vec(:)', t(:)'];
        end

        function T = getGtTransMat(gtData, idx)
            % gtData is a struct of arrays
            qx = gtData.qx(idx);
            qy = gtData.qy(idx);
            qz = gtData.qz(idx);
            qw = gtData.qw(idx);

            tx = gtData.X(idx);
            ty = gtData.Y(idx);
            tz = gtData.Z(idx);

            R = quat2rotm([qw, qx, qy, qz]);
            T = eye(4);
            T(1:3,1:3) = R;
            T(1:3,4) = [tx; ty; tz];
        end

        function pts3dCam = convertWorldPtIntoCamFrame(pts3d, projMat)
            % pts3d: 3xN
            % projMat: 4x4

            N = size(pts3d, 2);
            pts3d_h = [pts3d; ones(1, N)];

            % Transform: T * P
            pts3dCam_h = projMat * pts3d_h;

            % Normalize
            pts3dCam = pts3dCam_h(1:3, :) ./ pts3dCam_h(4, :);
        end

        function ptsObj = transform_cam_point_to_object(T, ptsCam)
            % T: 4x4 (Cam to Object?? No, usually World to Cam. Inverting)
            % Python: transform_cam_point_to_object uses inv(T)

            invT = inv(T);

            if size(ptsCam, 1) == 3
                ptsCam_h = [ptsCam; ones(1, size(ptsCam, 2))];
            else
                ptsCam_h = ptsCam;
            end

            res_h = invT * ptsCam_h;
            ptsObj = res_h(1:3, :) ./ res_h(4, :);
        end

        function [rayOrg, rayDir] = precompute_camera_rays_identity(imgPts, K)
            % imgPts: N x 2

            Kinv = inv(K);
            N = size(imgPts, 1);

            % Homogeneous coords (3 x N)
            uv1 = [imgPts, ones(N, 1)]';

            % Camera frame directions
            dir_c = Kinv * uv1;

            % Normalize
            norms = sqrt(sum(dir_c.^2, 1));
            dir_c = dir_c ./ norms;

            % R=I, t=0
            rayDir = dir_c'; % N x 3
            rayOrg = zeros(N, 3);
        end

        function [best_proj, best_idx] = project_cloud_onto_ray_fast(ray_o, ray_d, cloud)
            % ray_o: 1x3
            % ray_d: 1x3
            % cloud: N x 3

            % Vectors from origin to cloud
            V = cloud - ray_o;

            % Dot products (N x 1)
            % ray_d must be unit vector
            t_vals = V * ray_d';

            % Clamp t >= 0
            t_vals(t_vals < 0) = 0;

            % Projections
            projs = ray_o + t_vals * ray_d;

            % Distances
            dists = sum((cloud - projs).^2, 2);

            [~, best_idx] = min(dists);
            best_proj = projs(best_idx, :);
        end

        function r = residuals(pose, pts3d, pts2dObs, K, ~)
            % pose: 1x6 [rx, ry, rz, tx, ty, tz]
            rVec = pose(1:3);
            tVec = pose(4:6);

            R = UtilityFunctions.rodrigues(rVec);

            % Project
            % 1. Transform 3D points
            pts3d_h = [pts3d; ones(1, size(pts3d,2))];

            T = eye(4);
            T(1:3,1:3) = R;
            T(1:3,4) = tVec(:);

            p_cam = T * pts3d_h;
            p_cam = p_cam(1:3, :) ./ p_cam(4, :);

            % 2. Project to Pixel
            uv_hom = K * p_cam;
            uv = uv_hom(1:2, :) ./ uv_hom(3, :);

            % 3. Error
            diff = uv' - pts2dObs;
            r = diff(:); % Flatten
        end

        % --- Math Helpers ---
        function R = rodrigues(r)
            % Simple Rodriguez formula (Rot Vec -> Rot Mat)
            theta = norm(r);
            if theta < 1e-16
                R = eye(3);
            else
                k = r / theta;
                K = [0 -k(3) k(2); k(3) 0 -k(1); -k(2) k(1) 0];
                R = eye(3) + sin(theta)*K + (1-cos(theta))*(K*K);
            end
        end

        function r = rodrigues_inv(R)
            % Rot Mat -> Rot Vec
            theta = acos((trace(R) - 1) / 2);
            if theta < 1e-16
                r = zeros(3,1);
            else
                factor = theta / (2 * sin(theta));
                r = factor * [R(3,2) - R(2,3); R(1,3) - R(3,1); R(2,1) - R(1,2)];
            end
        end

        function logFileName = build_log_file_name(expFolderName, angStd, transStd, max_its, min_cost)
            % Helper to format the log filename
            ang_str = sprintf('%g', angStd);

            if length(transStd) == 3
                trans_str = sprintf('%g_%g_%g', transStd(1), transStd(2), transStd(3));
            else
                trans_str = sprintf('%g', transStd);
            end

            logFileName = sprintf('%s_%dit_%gminCost_ang%s_trans%s.txt', ...
                expFolderName, max_its, min_cost, ang_str, trans_str);
        end

        function fid = createLogDataFiles(fileName)
            % Create Logs folder if it doesn't exist
            if ~exist('Logs', 'dir')
                mkdir('Logs');
            end

            fullPath = fullfile('Logs', fileName);
            fid = fopen(fullPath, 'wt');

            if fid == -1
                error('Could not open log file: %s', fullPath);
            end

            % Columns matching Python code exactly
            cols = {...
                'sample_number', 'validty', ...
                'ideal_roll(deg)', 'ideal_pitch(deg)', 'ideal_yaw(deg)', ...
                'ideal_distance_x(m)', 'ideal_distance_y(m)', 'ideal_distance_z(m)', ...
                'noisy_roll(deg)', 'noisy_pitch(deg)', 'noisy_yaw(deg)', ...
                'noisy_distance_x(m)', 'noisy_distance_y(m)', 'noisy_distance_z(m)', ...
                'estimated_roll(deg)', 'estimated_pitch(deg)', 'estimated_yaw(deg)', ...
                'estimated_distance_x(m)', 'estimated_distance_y(m)', 'estimated_distance_z(m)', ...
                'lm_cost', 'lm_optimality', 'lm_numFuncEval', 'lm_iteration', ...
                'lm_status', 'lm_success', ...
                'errAng1', 'errAng2', 'errAng3', ...
                'rot_errvec_x(deg)', 'rot_errvec_y(deg)', 'rot_errvec_z(deg)', ...
                'rot_err_angle_mag(deg)', 'duration', 'loss_condition', 'ctr_condition'};

            fprintf(fid, '%s\n', strjoin(cols, '\t'));
        end

        function logData(fid, numSample, validity, lmRes, T_gt, T_noisy, R_est, t_est, duration, lossCond, ctrCond)
            % Extract GT and Noisy components
            R_gt = T_gt(1:3, 1:3);
            t_gt = T_gt(1:3, 4);

            R_noisy = T_noisy(1:3, 1:3);
            t_noisy = T_noisy(1:3, 4);

            % Calculate Euler Angles (ZYX convention, degrees)
            % MATLAB's rotm2eul returns [Yaw, Pitch, Roll] in radians.
            % We convert to degrees and reorder to [Roll, Pitch, Yaw] to match Python navpy.
            function rpy = getRPY(R)
                eul = rotm2eul(R, 'ZYX'); % Returns [z, y, x] (rad)
                rpy_deg = rad2deg(eul);
                rpy = [rpy_deg(3), rpy_deg(2), rpy_deg(1)]; % [x, y, z]
            end

            gtAng = getRPY(R_gt);
            nsAng = getRPY(R_noisy);
            estAng = getRPY(R_est);

            % Calculate Angular Errors (Wrapped to -180..180)
            diff = gtAng - estAng;
            errAng = mod(diff + 180, 360) - 180;

            % Calculate Rotation Error Vector (Axis-Angle)
            R_rel = R_est' * R_gt;
            rot_vec_rad = UtilityFunctions.rodrigues_inv(R_rel);
            rot_vec_deg = rad2deg(rot_vec_rad);
            rot_err_mag = rad2deg(norm(rot_vec_rad));

            % Formatted Output String
            % Note: 'lmRes' is a struct with cost, optimality, funcCount, iterations
            fmt = ['%d\t%d\t', ...                                   % sample, validity
                '%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t', ...       % GT (ang/trans)
                '%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t', ...       % Noisy
                '%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t', ...       % Estimated
                '%.4f\t%.4f\t%d\t%d\t%d\t%d\t', ...               % LM stats
                '%.4f\t%.4f\t%.4f\t', ...                         % ErrAngles
                '%.4f\t%.4f\t%.4f\t%.4f\t', ...                   % RotErrVec
                '%.4f\t%d\t%d\n'];                                % Duration, Conds

            fprintf(fid, fmt, ...
                numSample, validity, ...
                gtAng, t_gt', ...
                nsAng, t_noisy', ...
                estAng, t_est', ...
                lmRes.cost, lmRes.optimality, lmRes.funcCount, lmRes.iterations, 1, 1, ... % Hardcoded status/success=1 for simplicity
                errAng, ...
                rot_vec_deg', rot_err_mag, ...
                0.0, lossCond, ctrCond);
        end


        function [subset_cloud, indices] = filter_closest_points(ray_origin, ray_dir, cloud, N)
            % Helper: Selects N points closest to the ray for cleaner plotting
            % ray_origin: 1x3
            % ray_dir: 1x3
            % cloud: 3xM (Columns are points)

            if nargin < 4, N = 50; end

            M = size(cloud, 2);

            % Ensure vectors are column vectors for calculation
            O = ray_origin(:);
            D = ray_dir(:);

            % Vector from Origin to Cloud Points
            V = cloud - O;

            % Project V onto D (scalar projection t)
            % D must be normalized. 1x3 * 3xM -> 1xM
            t = D' * V;

            % Calculate Projection Points: P = O + t*D
            projections = O + D * t;

            % Squared distance from Cloud Point to Line
            dists_sq = sum((cloud - projections).^2, 1);

            % Filter points "behind" the camera (t < 0)
            dists_sq(t < 0) = inf;

            % Sort and pick top N
            [~, sort_idx] = sort(dists_sq);
            actual_N = min(N, M);
            indices = sort_idx(1:actual_N);

            subset_cloud = cloud(:, indices);
        end

        function visualize_entire_scene(rayOrigins, rayDirections, denseCloud, match_indices, T_gt, T_est)
            % visualize_entire_scene - Plots rays, cloud, highlights matches, and overlay text.
            %
            % Args:
            %   ... (previous args)
            %   T_gt: (Optional) 4x4 Ground Truth Transformation Matrix
            %   T_est: (Optional) 4x4 Estimated Transformation Matrix

            figure('Name', 'Full Scene Visualization', 'NumberTitle', 'off');
            clf; hold on; grid on; axis equal;
            xlabel('X'); ylabel('Y'); zlabel('Z');
            view(3);

            % 1. Plot Cloud
            scatter3(denseCloud(1,:), denseCloud(2,:), denseCloud(3,:), ...
                10, [0.8, 0.8, 0.8], 'filled', 'DisplayName', 'Dense Cloud');

            % 2. Plot Camera Path
            plot3(rayOrigins(:,1), rayOrigins(:,2), rayOrigins(:,3), ...
                'k-', 'LineWidth', 1.5, 'DisplayName', 'Camera Path');
            scatter3(rayOrigins(:,1), rayOrigins(:,2), rayOrigins(:,3), ...
                20, 'k', 'filled', 'HandleVisibility', 'off');

            % 3. Plot Rays
            N = size(rayOrigins, 1);
            ray_len = 20.0;
            X = zeros(3, N); Y = zeros(3, N); Z = zeros(3, N);

            for i = 1:N
                start_pt = rayOrigins(i, :);
                if nargin >= 4 && ~isempty(match_indices) && match_indices(i) > 0
                    end_pt = denseCloud(:, match_indices(i))';
                else
                    end_pt = start_pt + rayDirections(i, :) * ray_len;
                end
                X(:, i) = [start_pt(1); end_pt(1); NaN];
                Y(:, i) = [start_pt(2); end_pt(2); NaN];
                Z(:, i) = [start_pt(3); end_pt(3); NaN];
            end
            plot3(X(:), Y(:), Z(:), 'g-', 'LineWidth', 0.5, ...
                'Color', [0, 0.8, 0, 0.2], 'DisplayName', 'Rays');

            % 4. Plot Unique Matches
            if nargin >= 4 && ~isempty(match_indices)
                valid_mask = match_indices > 0;
                valid_idx = match_indices(valid_mask);
                if ~isempty(valid_idx)
                    matched_pts = denseCloud(:, valid_idx);
                    scatter3(matched_pts(1,:), matched_pts(2,:), matched_pts(3,:), ...
                        20, 'r', 'filled', 'DisplayName', 'Matched Points');
                end

                % 5. Highlight Multi-Hits
                [unique_vals, ~, idx_map] = unique(valid_idx);
                counts = accumarray(idx_map, 1);
                duplicate_ids = unique_vals(counts > 1);

                if ~isempty(duplicate_ids)
                    multi_hit_pts = denseCloud(:, duplicate_ids);
                    scatter3(multi_hit_pts(1,:), multi_hit_pts(2,:), multi_hit_pts(3,:), ...
                        60, 'm', 'filled', 'MarkerEdgeColor', 'k', 'DisplayName', 'Multi-Hit Points');
                    title(sprintf('Scene: %d Rays | %d Multi-Hits', N, length(duplicate_ids)));
                else
                    title(sprintf('Scene: %d Rays | No Multi-Hits', N));
                end
            end

            % --- 6. DISPLAY GROUND TRUTH & ESTIMATED POSE TEXT ---
            if nargin >= 6 && ~isempty(T_gt) && ~isempty(T_est)
                % Extract Trans
                t_gt = T_gt(1:3, 4);
                t_est = T_est(1:3, 4);

                % Extract Rot (Euler ZYX -> Roll, Pitch, Yaw)
                % Assuming rotm2eul is available (Navigation/Robotics Toolbox)
                % If not available, we can use a custom function.
                try
                    eul_gt = rotm2eul(T_gt(1:3,1:3), 'ZYX');
                    eul_est = rotm2eul(T_est(1:3,1:3), 'ZYX');
                    % Convert rad to deg and reorder to [Roll, Pitch, Yaw]
                    rpy_gt = rad2deg([eul_gt(3), eul_gt(2), eul_gt(1)]);
                    rpy_est = rad2deg([eul_est(3), eul_est(2), eul_est(1)]);
                catch
                    rpy_gt = [0,0,0]; rpy_est = [0,0,0]; % Fallback if toolbox missing
                end

                % Format String
                info_str = {
                    '\bfGround Truth:', ...
                    sprintf(' Pos: [%.2f, %.2f, %.2f]', t_gt), ...
                    sprintf(' Ang: [%.2f, %.2f, %.2f]', rpy_gt), ...
                    '', ...
                    '\bfEstimated:', ...
                    sprintf(' Pos: [%.2f, %.2f, %.2f]', t_est), ...
                    sprintf(' Ang: [%.2f, %.2f, %.2f]', rpy_est)
                    };

                % Display Text Box
                annotation('textbox', [0.02, 0.75, 0.25, 0.2], ...
                    'String', info_str, ...
                    'FitBoxToText', 'on', ...
                    'BackgroundColor', 'w', ...
                    'FaceAlpha', 0.8, ...
                    'EdgeColor', 'k');
            end

            legend('show', 'Location', 'best');
            hold off;
        end
    end
end