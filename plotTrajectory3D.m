function plotTrajectory3D(fileName)
    % plotTrajectory3D - Visualizes GT, Noisy, and Estimated paths in 3D.
    % Usage: plotTrajectory3D('path/to/logfile.txt')

    % 1. Read the Log File
    opts = detectImportOptions(fileName, 'FileType', 'text', 'Delimiter', '\t', 'NumHeaderLines', 1);
    data = readmatrix(fileName, opts);

    % 2. Extract Coordinates based on your log structure
    % Ground Truth (Cols 6,7,8)
    gtX = data(:, 6);
    gtY = data(:, 7);
    gtZ = data(:, 8);

    % Noisy Input (Cols 12,13,14)
    noisyX = data(:, 12);
    noisyY = data(:, 13);
    noisyZ = data(:, 14);

    % Estimated Output (Cols 18,19,20)
    estX = data(:, 18);
    estY = data(:, 19);
    estZ = data(:, 20);

    % 3. Create the Plot
    figure('Name', '3D Trajectory Reconstruction', 'Color', 'w', 'NumberTitle', 'off');
    hold on; grid on; axis equal;
    
    % --- Plot A: Ground Truth (Green, Solid) ---
    plot3(gtX, gtY, gtZ, 'g-', 'LineWidth', 2, 'DisplayName', 'Ideal (GT)');
    
    % --- Plot B: Noisy Input (Blue, Dotted, Semi-transparent) ---
    % We use a thinner line or dots so it doesn't obscure the result
    plot3(noisyX, noisyY, noisyZ, 'b:', 'LineWidth', 1, 'DisplayName', 'Noisy Input');
    
    % --- Plot C: Estimated (Red, Dashed) ---
    plot3(estX, estY, estZ, 'r--', 'LineWidth', 2, 'DisplayName', 'Estimated');

    % 4. Add Visual Aids
    % Mark Start Point
    scatter3(gtX(1), gtY(1), gtZ(1), 100, 'g', 'filled', 'MarkerEdgeColor', 'k', 'DisplayName', 'Start Point');
    
    % Optional: Draw "Error Lines" connecting Noisy -> Estimated for a few points
    % This visualizes the "correction" your algorithm performed.
    step = max(1, floor(length(gtX)/50)); % Draw ~50 lines total to avoid clutter
    for k = 1:step:length(gtX)
        plot3([noisyX(k), estX(k)], [noisyY(k), estY(k)], [noisyZ(k), estZ(k)], ...
              'Color', [0.7 0.7 0.7 0.5], 'LineWidth', 0.5, 'HandleVisibility', 'off');
    end

    % 5. Formatting
    xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
    title('Trajectory Reconstruction: Ideal vs Noisy vs Estimated');
    legend('show', 'Location', 'best');
    view(3); % Set default 3D view
    
    hold off;
end