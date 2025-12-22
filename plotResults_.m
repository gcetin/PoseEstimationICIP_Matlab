function stats = plotResults(fName, plot_all_angles)
    % plotResults - Visualizes Pose Estimation Logs (MATLAB Port)
    %
    % Args:
    %   fName (string): Path to the log file.
    %   plot_all_angles (bool): If true, plots Phi/Psi in addition to Theta.
    
    if nargin < 2
        plot_all_angles = false;
    end
    
    % close all;
    
    % --- 1. Filename Parsing ---
    [~, name, ~] = fileparts(fName);
    declared_samples = [];
    
    % Regex to find sample count (e.g., _500_IDEAL)
    tokens = regexp(name, '_(\d+)_(IDEAL|NOISY)', 'tokens');
    if ~isempty(tokens)
        declared_samples = str2double(tokens{1}{1});
        fprintf('Samples declared in filename: %d\n', declared_samples);
    end

    % --- 2. Load Data ---
    % 1 header line, tab delimiter
    opts = detectImportOptions(fName, 'FileType', 'text', 'Delimiter', '\t', 'NumHeaderLines', 1);
    raw_data = readmatrix(fName, opts);
    
    [total_rows, num_cols] = size(raw_data);
    raw_sampleNo = raw_data(:, 1);

    % --- 3. Averaging Logic ---
    block_size = length(raw_sampleNo);
    if ~isempty(declared_samples)
        block_size = declared_samples;
    end
    
    use_average = ~isempty(declared_samples) && mod(total_rows, block_size) == 0;
    repeats_info = 'unknown';

    if use_average
        repeat_cnt = total_rows / block_size;
        
        % Reshape data to [block_size, repeat_cnt, num_cols] to average across repeats
        % We transpose raw_data first to handle linear indexing correctly during reshape
        % Reshape logic: (Cols, BlockSize, Repeats) -> Permute to (BlockSize, Repeats, Cols)
        temp_data = reshape(raw_data', num_cols, block_size, repeat_cnt);
        perm_data = permute(temp_data, [2, 3, 1]);
        
        % Mean across the 2nd dimension (repeats)
        data = squeeze(mean(perm_data, 2));
        
        % Restore sample numbers (they shouldn't be averaged, but strictly they are identical)
        data(:, 1) = raw_sampleNo(1:block_size);
        
        repeats_info = repeat_cnt;
        fprintf('Averaged recordings per sample: block_size=%d, repeats=%d\n', block_size, repeat_cnt);
    else
        data = raw_data;
        fprintf('Averaging skipped: cannot infer repeat block.\n');
    end

    % --- 4. Column Extraction ---
    % MATLAB uses 1-based indexing, so Python index 0 is 1 here.
    sampleNo = data(:, 1);
    validty  = data(:, 2);

    % Ground Truth
    gtAngPhi   = data(:, 3);
    gtAngTheta = data(:, 4);
    gtAngPsi   = data(:, 5);
    gtX        = data(:, 6);
    gtY        = data(:, 7);
    gtZ        = data(:, 8);

    % Noisy Pose
    noisyAngPhi   = data(:, 9);
    noisyAngTheta = data(:, 10);
    noisyAngPsi   = data(:, 11);
    noisyX        = data(:, 12);
    noisyY        = data(:, 13);
    noisyZ        = data(:, 14);

    % Estimated Pose
    estAngPhi   = data(:, 15);
    estAngTheta = data(:, 16);
    estAngPsi   = data(:, 17);
    estX        = data(:, 18);
    estY        = data(:, 19);
    estZ        = data(:, 20);

    % LM Stats (Indices 21-26)
    lmCost      = data(:, 21);
    % ... others skipped as per python code ...

    % Errors (Indices 27-36)
    % errAng1 = data(:, 27); 
    rotErrMag      = data(:, 33);
    duration       = data(:, 34);
    loss_condition = data(:, 35);
    ctr_condition  = data(:, 36);

    % --- 5. Computations ---
    
    % Translation Differences
    err_est_x = estX - gtX;
    err_est_y = estY - gtY;
    err_est_z = estZ - gtZ;

    err_noisy_x = noisyX - gtX;
    err_noisy_y = noisyY - gtY;
    err_noisy_z = noisyZ - gtZ;

    % Percent Errors
    pe_est_x = percent_err(err_est_x, gtX);
    pe_est_y = percent_err(err_est_y, gtY);
    pe_est_z = percent_err(err_est_z, gtZ);
    pe_noisy_x = percent_err(err_noisy_x, gtX);
    pe_noisy_y = percent_err(err_noisy_y, gtY);
    pe_noisy_z = percent_err(err_noisy_z, gtZ);

    % Report Indices
    report_indices('X est-gt |error|>1', err_est_x, 1.0, sampleNo);
    report_indices('Y est-gt |error|>1', err_est_y, 1.0, sampleNo);
    report_indices('Z est-gt |error|>1', err_est_z, 5.0, sampleNo);

    % Angle Diffs
    phi_err   = angle_diff_deg(estAngPhi, gtAngPhi);
    theta_err = angle_diff_deg(estAngTheta, gtAngTheta);
    psi_err   = angle_diff_deg(estAngPsi, gtAngPsi);

    % % Angle Diffs
    % phi_err   = angle_diff_deg(abs(estAngPhi), abs(gtAngPhi));
    % theta_err = angle_diff_deg(abs(estAngTheta), abs(gtAngTheta));
    % psi_err   = angle_diff_deg(abs(estAngPsi), abs(gtAngPsi));    
    
    % --- 6. Plotting ---

    % % FIGURE 1: Translation Errors
    % fig1 = figure('Name', 'Translation Errors', 'NumberTitle', 'off');
    % set(gcf, 'Position', [100, 100, 1000, 800]);
    % 
    % ax1 = subplot(3,1,1); plot_trans_subplot(sampleNo, err_est_x, loss_condition, 'X-axis: translation differences');
    % ax2 = subplot(3,1,2); plot_trans_subplot(sampleNo, err_est_y, loss_condition, 'Y-axis: translation differences');
    % ax3 = subplot(3,1,3); plot_trans_subplot(sampleNo, err_est_z, loss_condition, 'Z-axis: translation differences');
    % xlabel('sample');
    % linkaxes([ax1, ax2, ax3], 'x');
    % 
    % % FIGURE 2: Percentage Translation Errors
    % fig2 = figure('Name', 'Percent Translation Errors', 'NumberTitle', 'off');
    % set(gcf, 'Position', [150, 150, 1000, 800]);
    % 
    % ax1 = subplot(3,1,1); plot_percent_subplot(sampleNo, pe_est_x, pe_noisy_x, 'X-axis: % translation error');
    % ax2 = subplot(3,1,2); plot_percent_subplot(sampleNo, pe_est_y, pe_noisy_y, 'Y-axis: % translation error');
    % ax3 = subplot(3,1,3); plot_percent_subplot(sampleNo, pe_est_z, pe_noisy_z, 'Z-axis: % translation error');
    % xlabel('sample');
    % linkaxes([ax1, ax2, ax3], 'x');

    % Setup Angle Data
    if plot_all_angles
        angle_names = {'Phi', 'Theta', 'Psi'};
        angle_errs = {phi_err, theta_err, psi_err};
        gt_angles = {gtAngPhi, gtAngTheta, gtAngPsi};
        est_angles = {estAngPhi, estAngTheta, estAngPsi};
        noisy_angles = {noisyAngPhi, noisyAngTheta, noisyAngPsi};
    else
        angle_names = {'Theta'};
        angle_errs = {theta_err};
        gt_angles = {gtAngTheta};
        est_angles = {estAngTheta};
        noisy_angles = {noisyAngTheta};
    end

    % FIGURE 3: Angle Errors
    figure('Name', 'Angle Errors (deg)', 'NumberTitle', 'off');
    num_ang = length(angle_names);
    axes_h = gobjects(num_ang, 1);
    for i = 1:num_ang
        axes_h(i) = subplot(num_ang, 1, i);
        plot(sampleNo, angle_errs{i}, '-o', 'LineWidth', 1.5, 'MarkerSize', 3);
        yline(0, 'k', 'LineWidth', 0.8);
        title([angle_names{i} ' Angle Error (deg)']);
        ylabel('Error (deg)');
        grid on;
    end
    xlabel('Frame');
    linkaxes(axes_h, 'x');

    % FIGURE 4: Angles GT vs Est
    figure('Name', 'Angles (GT vs Est)', 'NumberTitle', 'off');
    axes_h = gobjects(num_ang, 1);
    for i = 1:num_ang
        axes_h(i) = subplot(num_ang, 1, i);
        plot(sampleNo, gt_angles{i}, 'LineWidth', 1.1, 'DisplayName', 'gt'); hold on;
        plot(sampleNo, est_angles{i}, '--', 'LineWidth', 1.1, 'DisplayName', 'est');
        ylabel([angle_names{i} ' (deg)']);
        title([angle_names{i} ' angle: gt vs est']);
        grid on; legend show;
    end
    xlabel('Frame');
    linkaxes(axes_h, 'x');

    % FIGURE 5: Theta Error vs Z Translation
    figure('Name', 'Theta Error vs Z translation', 'NumberTitle', 'off');
    scatter(gtZ, theta_err, 12, 'filled', 'MarkerFaceAlpha', 0.85);
    xlabel('Ground-truth Z (m)');
    ylabel('Theta error (deg)');
    title('Theta error vs ground-truth Z translation');
    grid on;

    % FIGURE 6: Angle Errors vs |GT| (Sorted)
    figure('Name', 'Angle Errors vs |GT| (sorted)', 'NumberTitle', 'off');
    for i = 1:num_ang
        subplot(num_ang, 1, i);
        abs_gt = abs(gt_angles{i});
        [sorted_gt, idx] = sort(abs_gt);
        sorted_err = abs(angle_errs{i}(idx));
        plot(sorted_gt, sorted_err, '-o', 'LineWidth', 1.1, 'MarkerSize', 3);
        ylabel('|error| (deg)');
        title([angle_names{i} ': |error| vs |gt| (sorted)']);
        grid on;
    end
    xlabel('|gt| (deg)');

    % FIGURE 7: Angle Errors vs GT (Sorted, Signed)
    figure('Name', 'Angle Errors vs GT (sorted, signed)', 'NumberTitle', 'off');
    for i = 1:num_ang
        subplot(num_ang, 1, i);
        [sorted_gt, idx] = sort(gt_angles{i});
        sorted_err = angle_errs{i}(idx);
        plot(sorted_gt, sorted_err, '-o', 'LineWidth', 1.1, 'MarkerSize', 3);
        yline(0, 'k--', 'LineWidth', 0.8);
        ylabel('error (deg)');
        title([angle_names{i} ': error vs gt (sorted by gt, signed)']);
        grid on;
    end
    xlabel('gt (deg)');

    % FIGURE 8: Angles GT vs Noisy
    figure('Name', 'Angles (GT vs Noisy)', 'NumberTitle', 'off');
    axes_h = gobjects(num_ang, 1);
    for i = 1:num_ang
        axes_h(i) = subplot(num_ang, 1, i);
        plot(sampleNo, gt_angles{i}, 'LineWidth', 1.1, 'DisplayName', 'gt'); hold on;
        plot(sampleNo, noisy_angles{i}, '--', 'LineWidth', 1.1, 'DisplayName', 'noisy');
        ylabel([angle_names{i} ' (deg)']);
        title([angle_names{i} ' angle: gt vs noisy']);
        grid on; legend show;
    end
    xlabel('Frame');
    linkaxes(axes_h, 'x');

    % FIGURE 9: Translation Values (Sorted by GT)
    figure('Name', 'Translation Values (sorted by GT)', 'NumberTitle', 'off');
    ax1 = subplot(3,1,1); plot_trans_values_sorted(gtX, estX, noisyX, 'X');
    ax2 = subplot(3,1,2); plot_trans_values_sorted(gtY, estY, noisyY, 'Y');
    ax3 = subplot(3,1,3); plot_trans_values_sorted(gtZ, estZ, noisyZ, 'Z');
    xlabel('gt translation (m)');

    % FIGURE 10: Absolute Translation Errors vs Sorted GT
    figure('Name', 'Absolute Translation Errors vs Sorted GT', 'NumberTitle', 'off');
    ax1 = subplot(3,1,1); plot_trans_errors_sorted(gtX, err_est_x, 'X');
    ax2 = subplot(3,1,2); plot_trans_errors_sorted(gtY, err_est_y, 'Y');
    ax3 = subplot(3,1,3); plot_trans_errors_sorted(gtZ, err_est_z, 'Z');
    xlabel('gt translation (m)');

    % --- 7. Statistics & Console Output ---
    
    % Helper for stats
    get_stats = @(x, abs_flag) struct('min', min(x, [], 'omitnan'), ...
                                      'max', max(x, [], 'omitnan'), ...
                                      'mean', mean(x, 'omitnan'), ...
                                      'median', median(x, 'omitnan'));
    
    calc_print_stats = @(arr, name, is_abs) print_stats_helper(arr, name, is_abs, get_stats);

    % Translation Absolute
    fprintf('\n=== Translation Error (absolute, meters) ===\n');
    stats_trans_abs.x = calc_print_stats(abs(err_est_x), 'X', false);
    stats_trans_abs.y = calc_print_stats(abs(err_est_y), 'Y', false);
    stats_trans_abs.z = calc_print_stats(abs(err_est_z), 'Z', false);
    trans_err_mag = sqrt(err_est_x.^2 + err_est_y.^2 + err_est_z.^2);
    stats_trans_abs.mag = calc_print_stats(trans_err_mag, '3D magnitude', false);

    % Angle Absolute
    fprintf('\n=== Angle Error (absolute, degrees) ===\n');
    stats_ang_abs.phi = calc_print_stats(abs(phi_err), 'Phi', false);
    stats_ang_abs.theta = calc_print_stats(abs(theta_err), 'Theta', false);
    stats_ang_abs.psi = calc_print_stats(abs(psi_err), 'Psi', false);

    % Translation Signed
    fprintf('\n=== Translation Error (signed, meters; stats on abs) ===\n');
    stats_trans_signed.x = calc_print_stats(err_est_x, 'X', true);
    stats_trans_signed.y = calc_print_stats(err_est_y, 'Y', true);
    stats_trans_signed.z = calc_print_stats(err_est_z, 'Z', true);

    % Angle Signed
    fprintf('\n=== Angle Error (signed, degrees; stats on abs) ===\n');
    stats_ang_signed.phi = calc_print_stats(phi_err, 'Phi', true);
    stats_ang_signed.theta = calc_print_stats(theta_err, 'Theta', true);
    stats_ang_signed.psi = calc_print_stats(psi_err, 'Psi', true);

    % Rot Mag
    fprintf('\n=== Axis-angle Rotation Error (deg) ===\n');
    stats_rot.mag = calc_print_stats(rotErrMag, 'Magnitude', true);

    % Build Output Struct
    stats.samples_declared = declared_samples;
    stats.samples_loaded = total_rows;
    stats.samples_after_average = size(data, 1);
    stats.translation_abs = stats_trans_abs;
    stats.angle_abs = stats_ang_abs;
    stats.translation_signed = stats_trans_signed;
    stats.angle_signed = stats_ang_signed;

end

% --- Helper Functions ---

function err = percent_err(err, base)
    denom = max(abs(base), 1e-6);
    err = 100.0 * abs(err) ./ denom;
end

function d = angle_diff_deg(est, gt)
    % Normalize angle differences to [-180, 180)
    d = est - gt;
    d = mod(d + 180, 360) - 180;
end

function report_indices(label, err, threshold, sampleNo)
    mask = abs(err) > threshold;
    idx = sampleNo(mask);
    fprintf('%s: count=%d\n', label, length(idx));
    if ~isempty(idx)
        fprintf('  indices (sampleNo): %s\n', mat2str(idx', 5)); % mat2str handles listing
    end
end

function plot_trans_subplot(x, y, loss, titlestr)
    yyaxis left;
    plot(x, y, 'LineWidth', 1.0, 'DisplayName', 'est - gt');
    yline(0, 'k', 'LineWidth', 0.8);
    ylabel('error (m)');
    title(titlestr);
    grid on;
    
    yyaxis right;
    plot(x, loss, 'Color', [0.85 0.325 0.098], 'LineWidth', 0.9, 'DisplayName', 'loss cond');
    ylabel('loss condition');
    ylim([-0.1, 1.1]);
    
    % Set color of right axis to match plot
    ax = gca;
    ax.YAxis(2).Color = [0.85 0.325 0.098];
end

function plot_percent_subplot(x, y_est, y_noisy, titlestr)
    plot(x, y_est, 'Color', [0.4940 0.1840 0.5560], 'LineWidth', 0.9, 'DisplayName', '% error est'); hold on;
    plot(x, y_noisy, 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 0.9, 'DisplayName', '% error noisy');
    grid on;
    ylabel('% error');
    title(titlestr);
    ylim([0, 200]);
    legend show;
end

function plot_trans_values_sorted(gt, est, noisy, axis_name)
    [gt_s, idx] = sort(gt);
    est_s = est(idx);
    noisy_s = noisy(idx);
    
    plot(gt_s, est_s, 'DisplayName', 'est vs gt', 'LineWidth', 1.0); hold on;
    plot(gt_s, noisy_s, 'DisplayName', 'noisy vs gt', 'LineWidth', 1.0);
    plot([min(gt_s) max(gt_s)], [min(gt_s) max(gt_s)], 'k--', 'LineWidth', 0.9, 'DisplayName', 'y=x (gt)');
    
    grid on;
    ylabel('translation (m)');
    title([axis_name '-axis: translation values (sorted by gt)']);
    legend('Location', 'best');
end

function plot_trans_errors_sorted(gt, err_est, axis_name)
    [gt_s, idx] = sort(gt);
    err_est_s = abs(err_est(idx));
    
    plot(gt_s, err_est_s, 'DisplayName', 'abs(est - gt)', 'LineWidth', 1.0);
    yline(0, 'k--', 'LineWidth', 0.9);
    grid on;
    ylabel('abs error (m)');
    title([axis_name '-axis: absolute errors vs sorted gt']);
    legend('Location', 'best');
end

function s = print_stats_helper(arr, name, use_abs_for_calc, stats_fn)
    % Calculates stats, prints them formatted, and returns struct
    val_for_calc = arr;
    if use_abs_for_calc
        val_for_calc = abs(arr);
    end
    s = stats_fn(val_for_calc, false);
    fprintf('%s: min=%.4f, max=%.4f, mean=%.4f, median=%.4f\n', ...
        name, s.min, s.max, s.mean, s.median);
end