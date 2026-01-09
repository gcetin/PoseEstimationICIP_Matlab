function yaw_ambiguity_reprojection_demo()
    clc; close all;

    % -----------------------------
    % 1) Camera intrinsics K
    % -----------------------------
    % TODO: Replace this with *your* real K
    fx = 1500; fy = 1500; cx = 960; cy = 540;  % example for 1920x1080
    K = [fx 0  cx;
         0  fy cy;
         0  0  1];

    imgW = 1920; imgH = 1080;

    % -----------------------------
    % 2) Generate 3D eight-shape points
    %    (replace with your exact model if you want)
    % -----------------------------
    N = 400;
    t = linspace(0, 2*pi, N);

    % Classic lemniscate-ish param (in meters)
    % Feel free to adjust A,B to match your model scale.
    A = 10;   % lateral amplitude (m)
    B = 6;    % vertical amplitude (m)

    % 3D curve in target frame (Z=0 plane here; you can make it non-planar if you want)
    X = A * sin(t);
    Y = B * sin(t).*cos(t);
    Z = zeros(size(t));

    P = [X(:) Y(:) Z(:)];  % Nx3

    % -----------------------------
    % 3) Pose setup
    % -----------------------------
    z_translation = 100;                % meters
    tvec = [0; 0; z_translation];       % x=y=0, z=100m
    theta_deg = 45;                     % choose test yaw
    theta = deg2rad(theta_deg);

    R_pos = angle2dcm(0, theta, 0);
    R_neg = angle2dcm(0, -theta, 0);

    % R_pos = Rz(theta);   % yaw +theta
    % R_neg = Rz(-theta);  % yaw -theta

    % Transform to camera frame: Pc = R * P + t
    Pc_pos = (R_pos * P') + tvec;
    Pc_neg = (R_neg * P') + tvec;

    % -----------------------------
    % 4) Project
    % -----------------------------
    uv_pos = project_points(K, Pc_pos); % Nx2
    uv_neg = project_points(K, Pc_neg); % Nx2

    % Optional: keep only points in front of camera (Z>0)
    valid = (Pc_pos(3,:)' > 1e-6) & (Pc_neg(3,:)' > 1e-6);
    uv_pos = uv_pos(valid,:);
    uv_neg = uv_neg(valid,:);

    % Optional: keep only points inside image
    in_pos = uv_pos(:,1)>=1 & uv_pos(:,1)<=imgW & uv_pos(:,2)>=1 & uv_pos(:,2)<=imgH;
    in_neg = uv_neg(:,1)>=1 & uv_neg(:,1)<=imgW & uv_neg(:,2)>=1 & uv_neg(:,2)<=imgH;
    in = in_pos & in_neg;
    uv_pos = uv_pos(in,:);
    uv_neg = uv_neg(in,:);

    % -----------------------------
    % 5) Compare reprojections
    % -----------------------------
    d = sqrt(sum((uv_pos - uv_neg).^2, 2));   % per-point pixel distance

    rmse = sqrt(mean(d.^2));
    mean_d = mean(d);
    max_d = max(d);

    fprintf('Yaw +%.1f vs -%.1f reprojection difference:\n', theta_deg, theta_deg);
    fprintf('  Mean |Δu,v|: %.4f px\n', mean_d);
    fprintf('  RMSE |Δu,v|: %.4f px\n', rmse);
    fprintf('  Max  |Δu,v|: %.4f px\n', max_d);
    fprintf('  Points compared: %d\n', size(uv_pos,1));

    % -----------------------------
    % 6) Plots
    % -----------------------------
    figure('Name','Overlay of Projections'); hold on; grid on;
    plot(uv_pos(:,1), uv_pos(:,2), '.', 'MarkerSize', 10);
    plot(uv_neg(:,1), uv_neg(:,2), '.', 'MarkerSize', 10);
    set(gca,'YDir','reverse');
    xlabel('u (px)'); ylabel('v (px)');
    title(sprintf('Projection overlay: yaw +%.1f° vs -%.1f°', theta_deg, theta_deg));
    legend('+yaw','-yaw');

    figure('Name','Per-point Pixel Difference'); grid on;
    plot(d, 'LineWidth', 1.5);
    xlabel('Point index along trajectory');
    ylabel('|Δ(u,v)| (px)');
    title(sprintf('Per-point reprojection difference (RMSE=%.4f px)', rmse));

    figure('Name','Histogram of Pixel Differences'); grid on;
    histogram(d, 40);
    xlabel('|Δ(u,v)| (px)');
    ylabel('Count');
    title('Distribution of yaw-sign reprojection differences');
end

% --- helpers ---
function R = Rz(a)
    ca = cos(a); sa = sin(a);
    R = [ ca -sa 0;
          sa  ca 0;
          0   0  1 ];
end

function uv = project_points(K, Pc)
    % Pc: 3xN points in camera frame
    x = Pc(1,:)./Pc(3,:);
    y = Pc(2,:)./Pc(3,:);
    uvh = K * [x; y; ones(1, numel(x))];
    uv = uvh(1:2,:)';
end
