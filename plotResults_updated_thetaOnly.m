
function stats = plotResults_updated_thetaOnly(fName, plot_all_angles)
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

    % Angle Diff (theta only)
    theta_err = angle_diff_deg(estAngTheta, gtAngTheta);
% --- 6. Plotting (Meaningful Plot Pack) ---
    % Notes:
    %  - We prioritize distribution/tail behavior (CDF + percentiles),
    %    distance dependence (gtZ), and optimizer diagnostics (lmCost/runtime),
    %    plus noisy->estimated improvement.
    %
    % Mask invalid samples if validity flag exists.
    valid_mask = true(size(sampleNo));
    if exist('validty','var') && ~isempty(validty)
        % Many logs use 0 = valid. If this removes everything, fall back to all rows.
        valid_mask = valid_mask & (validty == 0);
        if ~any(valid_mask)
            warning('valid_mask empty; plotting all rows (check Validity convention).');
            valid_mask = true(size(sampleNo));
        end
    end

    % Translation error magnitudes
    e3D_est   = sqrt(err_est_x.^2 + err_est_y.^2 + err_est_z.^2);
    e3D_noisy = sqrt(err_noisy_x.^2 + err_noisy_y.^2 + err_noisy_z.^2);

    % ---------- Figure A: CDFs of translation errors ----------
    figA = figure('Name', 'CDF: Translation Errors', 'NumberTitle', 'off');
    tiledlayout(figA, 2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile; hold on; grid on;
    plot_cdf(abs(err_est_x(valid_mask)), '|ΔX| (est)');
    plot_cdf(abs(err_noisy_x(valid_mask)), '|ΔX| (noisy)');
    xlabel('meters'); ylabel('CDF'); title('X-axis'); legend('Location','southeast');

    nexttile; hold on; grid on;
    plot_cdf(abs(err_est_y(valid_mask)), '|ΔY| (est)');
    plot_cdf(abs(err_noisy_y(valid_mask)), '|ΔY| (noisy)');
    xlabel('meters'); ylabel('CDF'); title('Y-axis'); legend('Location','southeast');

    nexttile; hold on; grid on;
    plot_cdf(abs(err_est_z(valid_mask)), '|ΔZ| (est)');
    plot_cdf(abs(err_noisy_z(valid_mask)), '|ΔZ| (noisy)');
    xlabel('meters'); ylabel('CDF'); title('Z-axis'); legend('Location','southeast');

    nexttile; hold on; grid on;
    plot_cdf(e3D_est(valid_mask), '||ΔT|| (est)');
    plot_cdf(e3D_noisy(valid_mask), '||ΔT|| (noisy)');
    xlabel('meters'); ylabel('CDF'); title('3D Translation'); legend('Location','southeast');

    % ---------- Figure B: CDF of rotation errors ----------
    figB = figure('Name', 'CDF: Rotation Errors', 'NumberTitle', 'off');
    hold on; grid on;
    has_rotmag = exist('rotErrMag','var') && ~all(isnan(rotErrMag(valid_mask)));
    if has_rotmag
        plot_cdf(rotErrMag(valid_mask), 'rotErrMag (axis-angle)');
    end
    plot_cdf(abs(theta_err(valid_mask)), '|Δθ| (Euler)');
    xlabel('degrees'); ylabel('CDF');
    title('Rotation error distributions');
    legend('Location','southeast');

    % ---------- Figure C: Key percentiles for 3D translation ----------
    figC = figure('Name', 'Percentiles: 3D Translation Error', 'NumberTitle', 'off');
    p = [0.50 0.90 0.95 0.99];
    q = quantile_safe(e3D_est(valid_mask), p);
    bar(100*p, q);
    grid on;
    xlabel('percentile (%)'); ylabel('meters');
    title('3D translation error percentiles (estimated)');

    % ---------- Figure D: Histograms (3D translation + rotation magnitude) ----------
    figD = figure('Name', 'Histograms: Error Magnitudes', 'NumberTitle', 'off');
    tiledlayout(figD, 2, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile; histogram(e3D_est(valid_mask), 50);
    grid on; xlabel('||ΔT|| (m)'); ylabel('Count');
    title('3D translation error (estimated)');

    nexttile;
    if has_rotmag
        histogram(rotErrMag(valid_mask), 50);
        xlabel('rotErrMag (deg)');
        title('Rotation error magnitude (axis-angle)');
    else
        histogram(abs(theta_err(valid_mask)), 50);
        xlabel('|Δθ| (deg)');
        title('Rotation error (Euler θ)');
    end
    grid on; ylabel('Count');

    % ---------- Figure E: Error vs distance (gtZ) with binned median ----------
    if exist('gtZ','var') && ~all(isnan(gtZ(valid_mask)))
        figE = figure('Name', 'Error vs Distance (gtZ)', 'NumberTitle', 'off');
        tiledlayout(figE, 2, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

        nexttile; hold on; grid on;
        scatter(gtZ(valid_mask), e3D_est(valid_mask), 12, 'filled', 'MarkerFaceAlpha', 0.35);
        [bc, bm] = binned_median(gtZ(valid_mask), e3D_est(valid_mask), 12);
        plot(bc, bm, 'LineWidth', 2);
        xlabel('gtZ (m)'); ylabel('||ΔT|| (m)');
        title('3D translation error vs distance');

        nexttile; hold on; grid on;
        scatter(gtZ(valid_mask), err_est_z(valid_mask), 12, 'filled', 'MarkerFaceAlpha', 0.35);
        [bc2, bm2] = binned_median(gtZ(valid_mask), err_est_z(valid_mask), 12);
        plot(bc2, bm2, 'LineWidth', 2);
        xlabel('gtZ (m)'); ylabel('ΔZ (m) (signed)');
        title('Signed Z error vs distance');
    end

    % ---------- Figure F: Optimizer diagnostics (cost + runtime) ----------
    has_cost = exist('lmCost','var') && ~all(isnan(lmCost(valid_mask)));
    has_time = exist('duration','var') && ~all(isnan(duration(valid_mask)));
    if has_cost || has_time
        figF = figure('Name', 'Diagnostics: Cost/Runtime vs Error', 'NumberTitle', 'off');
        tiledlayout(figF, 2, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

        nexttile; hold on; grid on;
        if has_cost
            x = log10(lmCost(valid_mask) + 1e-12);
            scatter(x, e3D_est(valid_mask), 12, 'filled', 'MarkerFaceAlpha', 0.35);
            xlabel('log10(lmCost)'); ylabel('||ΔT|| (m)');
            title('Error vs LM cost');
        else
            text(0.1, 0.5, 'lmCost not available', 'Units','normalized');
            axis off;
        end

        nexttile; hold on; grid on;
        if has_time
            scatter(duration(valid_mask), e3D_est(valid_mask), 12, 'filled', 'MarkerFaceAlpha', 0.35);
            xlabel('duration (s)'); ylabel('||ΔT|| (m)');
            title('Error vs runtime');
        else
            text(0.1, 0.5, 'duration not available', 'Units','normalized');
            axis off;
        end
    end

    % ---------- Figure G: Noisy -> Estimated improvement ----------
    delta_3D = e3D_noisy - e3D_est; % positive => improvement
    figG = figure('Name', 'Improvement: Noisy -> Estimated', 'NumberTitle', 'off');
    tiledlayout(figG, 2, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile; histogram(delta_3D(valid_mask), 60);
    grid on; xlabel('\Delta 3D (m)'); ylabel('Count');
    title('Improvement distribution (positive means better than noisy)');

    nexttile; hold on; grid on;
    if exist('gtZ','var') && ~all(isnan(gtZ(valid_mask)))
        scatter(gtZ(valid_mask), delta_3D(valid_mask), 12, 'filled', 'MarkerFaceAlpha', 0.35);
        [bc3, bm3] = binned_median(gtZ(valid_mask), delta_3D(valid_mask), 12);
        plot(bc3, bm3, 'LineWidth', 2);
        xlabel('gtZ (m)'); ylabel('\Delta 3D (m)');
        title('Improvement vs distance');
    else
        plot(sampleNo(valid_mask), delta_3D(valid_mask), '-');
        xlabel('sample'); ylabel('\Delta 3D (m)');
        title('Improvement vs sample index');
    end


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
    fprintf('=== Theta Error (absolute, degrees) ===');
    stats_ang_abs.theta = calc_print_stats(abs(theta_err), 'Theta', false);

    % Translation Signed
    fprintf('\n=== Translation Error (signed, meters; stats on abs) ===\n');
    stats_trans_signed.x = calc_print_stats(err_est_x, 'X', true);
    stats_trans_signed.y = calc_print_stats(err_est_y, 'Y', true);
    stats_trans_signed.z = calc_print_stats(err_est_z, 'Z', true);

    % Angle Signed
    fprintf('=== Theta Error (signed, degrees; stats on abs) ===');
    stats_ang_signed.theta = calc_print_stats(theta_err, 'Theta', true);

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


function plot_cdf(x, labelStr)
    % plot_cdf - Plot empirical CDF of vector x (ignores NaNs)
    x = x(:);
    x = x(~isnan(x));
    if isempty(x)
        return;
    end
    x = sort(x);
    n = numel(x);
    y = (1:n) / n;
    plot(x, y, 'DisplayName', labelStr, 'LineWidth', 1.5);
end

function q = quantile_safe(x, p)
    % quantile_safe - Quantiles with NaN handling
    x = x(:);
    x = x(~isnan(x));
    if isempty(x)
        q = nan(size(p));
        return;
    end
    q = prctile(x, 100*p);
end

function [binCenters, binMedian] = binned_median(x, y, nbins)
    % binned_median - Median of y in x-bins (ignores NaNs)
    x = x(:); y = y(:);
    m = ~isnan(x) & ~isnan(y);
    x = x(m); y = y(m);
    if isempty(x)
        binCenters = []; binMedian = [];
        return;
    end
    xmin = min(x); xmax = max(x);
    if xmin == xmax
        binCenters = xmin;
        binMedian = median(y);
        return;
    end
    edges = linspace(xmin, xmax, nbins+1);
    binCenters = (edges(1:end-1) + edges(2:end)) / 2;
    binMedian = nan(1, nbins);
    for i = 1:nbins
        idx = x >= edges(i) & x < edges(i+1);
        if i == nbins
            idx = x >= edges(i) & x <= edges(i+1);
        end
        if any(idx)
            binMedian(i) = median(y(idx));
        end
    end
    keep = ~isnan(binMedian);
    binCenters = binCenters(keep);
    binMedian = binMedian(keep);
end
