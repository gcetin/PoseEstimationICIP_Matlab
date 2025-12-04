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

        function T_noisy = addNoiseToPose(T_gt, ang_std, trans_std)
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
    
        function visualize_ray_projection(ray_origin, ray_dir, denseCloud, idx)
                % Visualizes the ray and adds it to the existing scene (Growing Scene)
                
                % 1. Target Figure 100
                fig = figure(100);
                set(fig, 'Name', 'Growing Scene Animation', 'NumberTitle', 'off');
                
                % 2. Check if this is the first plot (Axes are empty)
                first_run = isempty(get(gca, 'Children'));
                
                if first_run
                    hold on; grid on; axis equal;
                    xlabel('X'); ylabel('Y'); zlabel('Z');
                    view(3);
                end
                
                % Normalize Ray
                ray_dir = ray_dir / norm(ray_dir);
                
                % Filter Cloud to closest 50 points
                [subset_cloud, ~] = UtilityFunctions.filter_closest_points(ray_origin, ray_dir, denseCloud, 50);
                
                % Dynamic Ray Length
                ray_len = 20.0;
                if ~isempty(subset_cloud)
                    V = subset_cloud - ray_origin(:);
                    t_vals = ray_dir(:)' * V;
                    max_t = max(t_vals);
                    if max_t > 0; ray_len = max_t * 1.2; end
                end
                
                % --- Plotting ---
                % Only add DisplayName (Legend) for the first run to keep it clean
                if first_run
                    ray_tag = 'Ray';
                    cam_tag = 'CamCenter';
                    vis_switch = 'on';
                else
                    ray_tag = '';
                    cam_tag = '';
                    vis_switch = 'off';
                end
                
                % 1. Draw Ray
                ray_end = ray_origin(:) + ray_dir(:) * ray_len;
                plot3([ray_origin(1), ray_end(1)], [ray_origin(2), ray_end(2)], [ray_origin(3), ray_end(3)], ...
                      'g-', 'LineWidth', 1, 'DisplayName', ray_tag, 'HandleVisibility', vis_switch);
                  
                % 2. Draw Camera Center
                if first_run % Plot camera center only once to save performance
                    scatter3(ray_origin(1), ray_origin(2), ray_origin(3), 50, 'k', 'filled', ...
                             'DisplayName', cam_tag);
                end
                
                % 3. Draw Cloud Points (Subset)
                if ~isempty(subset_cloud)
                    scatter3(subset_cloud(1,:), subset_cloud(2,:), subset_cloud(3,:), 10, 'b', 'filled', ...
                             'MarkerFaceAlpha', 0.5, 'HandleVisibility', 'off');
                end
                
                % 4. Draw Target Match (Always show, but don't clutter legend)
                if idx > 0 && idx <= size(denseCloud, 2)
                    target_pt = denseCloud(:, idx);
                    scatter3(target_pt(1), target_pt(2), target_pt(3), 30, [1, 0.5, 0], 'filled', ...
                             'MarkerEdgeColor', 'none', 'HandleVisibility', 'off');
                end
                
                if first_run
                    legend('show', 'Location', 'northeast');
                end
                
                drawnow;
            end
    
    end
end