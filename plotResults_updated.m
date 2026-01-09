
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

    % --- 3. Median Logic (Replaces Averaging) ---
    block_size = length(raw_sampleNo);
    if ~isempty(declared_samples)
        block_size = declared_samples;
    end
    
    use_average = ~isempty(declared_samples) && mod(total_rows, block_size) == 0;
    repeats_info = 'unknown';
    % use_average = false
    if use_average
        repeat_cnt = total_rows / block_size;
        
        % Reshape data to [block_size, repeat_cnt, num_cols]
        temp_data = reshape(raw_data', num_cols, block_size, repeat_cnt);
        perm_data = permute(temp_data, [2, 3, 1]);
        
        % --- CHANGE: Use Median instead of Mean ---
        % This filters out outliers across the repeats
        data = squeeze(median(perm_data, 2));
        
        % Restore sample numbers (strictly they are identical)
        data(:, 1) = raw_sampleNo(1:block_size);
        
        repeats_info = repeat_cnt;
        fprintf('Median-filtering recordings per sample: block_size=%d, repeats=%d\n', block_size, repeat_cnt);
    else
        data = raw_data;
        fprintf('Median logic skipped: cannot infer repeat block.\n');
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
    % --- 6. Plotting (Meaningful Plot Pack) ---
    %
    % Goal: fewer plots, more insight.
    %   A) Accuracy & tail risk: CDF + histogram
    %   B) Bias vs distance: error vs gtZ (or range)
    %   C) Convergence diagnostics: error vs LM cost / duration
    %   D) Improvement: noisy -> estimated

        % --- Derived metrics used by plots ---
        e_est_3D   = sqrt(err_est_x.^2 + err_est_y.^2 + err_est_z.^2);
        e_noisy_3D = sqrt(err_noisy_x.^2 + err_noisy_y.^2 + err_noisy_z.^2);

        abs_ex = abs(err_est_x); abs_ey = abs(err_est_y); abs_ez = abs(err_est_z);

        % Some logs may contain NaNs; filter those in the plots
        valid_mask = isfinite(sampleNo) & isfinite(e_est_3D) & isfinite(gtZ);

        % =========================
        % FIG 1: CDFs (accuracy + tails)
        % =========================
        fig1 = figure('Name','CDF: Translation & Rotation Errors','NumberTitle','off');
        set(fig1, 'Position', [100, 100, 1050, 750]);

        subplot(2,2,1);
        plot_cdf(abs_ex(valid_mask), '|\Delta X| (m)'); hold on;
        plot_cdf(abs_ey(valid_mask), '|\Delta Y| (m)');
        plot_cdf(abs_ez(valid_mask), '|\Delta Z| (m)');
        plot_cdf(e_est_3D(valid_mask), '3D error (m)');
        grid on; xlabel('Error'); ylabel('CDF'); title('Translation error CDF');
        legend('Location','southeast');

        subplot(2,2,2);
        if exist('rotErrMag','var') && ~isempty(rotErrMag)
            plot_cdf(abs(rotErrMag(valid_mask)), 'rotErrMag');
            grid on; xlabel('Rotation error'); ylabel('CDF'); title('Rotation error CDF (axis-angle)');
            legend('Location','southeast');
        else
            text(0.1,0.5,'rotErrMag not available','FontSize',12); axis off;
        end

        subplot(2,2,3);
        plot_cdf(abs(theta_err(valid_mask)), '|\Delta\theta| (deg)'); hold on;
        if plot_all_angles
            plot_cdf(abs(phi_err(valid_mask)),   '|\Delta\phi| (deg)');
            plot_cdf(abs(psi_err(valid_mask)),   '|\Delta\psi| (deg)');
        end
        grid on; xlabel('Angle error (deg)'); ylabel('CDF'); title('Euler angle error CDF');
        legend('Location','southeast');

        subplot(2,2,4);
        % Annotate key quantiles for 3D translation error
        qs = [0.50 0.90 0.95 0.99];
        qv = quantile_safe(e_est_3D(valid_mask), qs);
        bar(qs*100, qv); grid on;
        xlabel('Percentile'); ylabel('3D error (m)');
        title('3D translation error percentiles');

        % =========================
        % FIG 2: Histograms (distribution shape)
        % =========================
        fig2 = figure('Name','Hist: 3D Translation Error & rotErrMag','NumberTitle','off');
        set(fig2, 'Position', [160, 120, 1050, 450]);

        subplot(1,2,1);
        histogram(e_est_3D(valid_mask), 60); grid on;
        xlabel('3D translation error (m)'); ylabel('Count'); title('Distribution of 3D translation error');

        subplot(1,2,2);
        if exist('rotErrMag','var') && ~isempty(rotErrMag)
            histogram(abs(rotErrMag(valid_mask)), 60); grid on;
            xlabel('rotErrMag'); ylabel('Count'); title('Distribution of rotErrMag');
        else
            text(0.1,0.5,'rotErrMag not available','FontSize',12); axis off;
        end

        % =========================
        % FIG 3: Bias vs distance (gtZ) with binned median curve
        % =========================
        fig3 = figure('Name','Error vs gtZ (bias check)','NumberTitle','off');
        set(fig3, 'Position', [120, 140, 1050, 600]);

        subplot(2,1,1);
        scatter(gtZ(valid_mask), e_est_3D(valid_mask), 8, 'filled'); grid on;
        xlabel('GT Z (m)'); ylabel('3D translation error (m)'); title('3D error vs gtZ');
        hold on;
        [bx, bmed] = binned_median(gtZ(valid_mask), e_est_3D(valid_mask), 25);
        plot(bx, bmed, 'LineWidth', 2); % default color
        legend({'Samples','Binned median'}, 'Location','northwest');

        subplot(2,1,2);
        scatter(gtZ(valid_mask), err_est_z(valid_mask), 8, 'filled'); grid on;
        xlabel('GT Z (m)'); ylabel('\Delta Z = estZ - gtZ (m)'); title('Signed Z error vs gtZ');
        hold on;
        [bx2, bmed2] = binned_median(gtZ(valid_mask), err_est_z(valid_mask), 25);
        plot(bx2, bmed2, 'LineWidth', 2);
        yline(0, '--'); legend({'Samples','Binned median','0'}, 'Location','northwest');

        % =========================
        % FIG 4: Diagnostics (LM cost / duration)
        % =========================
        fig4 = figure('Name','Diagnostics: Error vs LM cost & duration','NumberTitle','off');
        set(fig4, 'Position', [140, 160, 1050, 600]);

        subplot(2,1,1);
        if exist('lmCost','var') && ~isempty(lmCost)
            x = log10(lmCost(valid_mask) + eps);
            scatter(x, e_est_3D(valid_mask), 8, 'filled'); grid on;
            xlabel('log10(lmCost)'); ylabel('3D translation error (m)');
            title('Error vs LM cost (should correlate if cost is meaningful)');
        else
            text(0.1,0.5,'lmCost not available','FontSize',12); axis off;
        end

        subplot(2,1,2);
        if exist('duration','var') && ~isempty(duration)
            scatter(duration(valid_mask), e_est_3D(valid_mask), 8, 'filled'); grid on;
            xlabel('Duration'); ylabel('3D translation error (m)');
            title('Error vs runtime (slow cases may indicate hard geometry)');
        else
            text(0.1,0.5,'duration not available','FontSize',12); axis off;
        end

        % =========================
        % FIG 5: Improvement (Noisy -> Estimated)
        % =========================
        fig5 = figure('Name','Improvement: Noisy -> Estimated (Translation)','NumberTitle','off');
        set(fig5, 'Position', [180, 180, 1050, 600]);

        delta_3D = e_noisy_3D - e_est_3D; % positive => improvement
        subplot(2,1,1);
        histogram(delta_3D(valid_mask), 60); grid on;
        xlabel('\Delta = noisy3D - est3D (m)'); ylabel('Count');
        title('Improvement distribution (positive means better than noisy)');

        subplot(2,1,2);
        scatter(gtZ(valid_mask), delta_3D(valid_mask), 8, 'filled'); grid on;
        xlabel('GT Z (m)'); ylabel('\Delta (m)'); title('Improvement vs gtZ');
        yline(0, '--');

        pct_improved = mean(delta_3D(valid_mask) > 0, 'omitnan') * 100;
        fprintf('\n[PlotPack] Translation improvement: %.2f%% of samples improved over noisy (delta_3D > 0).\n', pct_improved);

        % Console: worst-k samples by 3D error
        k = min(20, nnz(valid_mask));
        if k > 0
            sample_valid = sampleNo(valid_mask);
            gtZ_valid    = gtZ(valid_mask);
            e3d_valid    = e_est_3D(valid_mask);
            [~, sidx] = sort(e3d_valid, 'descend');
            fprintf('[PlotPack] Worst-%d samples by 3D translation error:\n', k);
            for ii = 1:k
                jj = sidx(ii);
                fprintf('  sample %d: e3D=%.4f m (gtZ=%.2f m)\n', sample_valid(jj), e3d_valid(jj), gtZ_valid(jj));
            end
        end

    xlabel('\Delta 3D (m)'); ylabel('Count');
    title('Improvement distribution (positive means better than noisy)');

    subplot(2,1,2);
    scatter(gtZ(valid_mask), delta_3D(valid_mask), 8, 'filled'); grid on;
    xlabel('GT Z (m)'); ylabel('\Delta (m)'); title('Improvement vs gtZ');
    yline(0, '--');

    pct_improved = mean(delta_3D(valid_mask) > 0, 'omitnan') * 100;
    fprintf('\n[PlotPack] Translation improvement: %.2f%% of samples improved over noisy (delta_3D > 0).\n', pct_improved);

    % Console: worst-k samples by 3D error
    k = min(20, nnz(valid_mask));
    if k > 0
        sample_valid = sampleNo(valid_mask);
        gtZ_valid    = gtZ(valid_mask);
        e3d_valid    = e_est_3D(valid_mask);
        [~, sidx] = sort(e3d_valid, 'descend');
        fprintf('[PlotPack] Worst-%d samples by 3D translation error:\n', k);
        for ii = 1:k
            jj = sidx(ii);
            fprintf('  sample %d: e3D=%.4f m (gtZ=%.2f m)\n', sample_valid(jj), e3d_valid(jj), gtZ_valid(jj));
        end
    end



    
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
    % Calculates stats, prints them formatted (Median focused), and returns struct
    val_for_calc = arr;
    if use_abs_for_calc
        val_for_calc = abs(arr);
    end
    s = stats_fn(val_for_calc, false);
    
    % --- CHANGE: Print Median first/prominently rather than Mean ---
    fprintf('%s: median=%.4f (mean=%.4f), min=%.4f, max=%.4f\n', ...
        name, s.median, s.mean, s.min, s.max);
end


function plot_cdf(x, label_str)
    % Plot empirical CDF of x (ignores NaNs/Infs)
    x = x(isfinite(x));
    if isempty(x)
        return;
    end
    x = sort(x(:));
    n = numel(x);
    y = (1:n) / n;
    plot(x, y, 'DisplayName', label_str, 'LineWidth', 1.4);
end

function qv = quantile_safe(x, qs)
    % Safe quantile wrapper (ignores NaNs/Infs, returns NaN if empty)
    x = x(isfinite(x));
    if isempty(x)
        qv = nan(size(qs));
        return;
    end
    try
        qv = quantile(x, qs);
    catch
        % Older MATLAB versions: use prctile
        qv = prctile(x, qs*100);
    end
end

function [bx, bmed] = binned_median(x, y, nbins)
    % Compute median(y) in nbins across x range. Returns bin centers bx.
    x = x(:); y = y(:);
    m = isfinite(x) & isfinite(y);
    x = x(m); y = y(m);
    if isempty(x)
        bx = []; bmed = [];
        return;
    end

    xmin = min(x); xmax = max(x);
    if xmin == xmax
        bx = xmin;
        bmed = median(y, 'omitnan');
        return;
    end

    edges = linspace(xmin, xmax, nbins+1);
    bx = 0.5 * (edges(1:end-1) + edges(2:end));
    bmed = nan(size(bx));

    for i = 1:nbins
        mask = (x >= edges(i)) & (x < edges(i+1));
        if any(mask)
            bmed(i) = median(y(mask), 'omitnan');
        end
    end

    % Drop empty bins
    keep = isfinite(bmed);
    bx = bx(keep);
    bmed = bmed(keep);
end
