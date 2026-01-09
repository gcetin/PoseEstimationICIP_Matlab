classdef TrajectoryGenerator
    methods (Static)
        function points = create_dense_trajectory(shapeType, appendFlag, dense_deltaT)
            SHAPE_SIZE = 5.0;
            t = 0:dense_deltaT:(2*pi);
            
            if shapeType == 0 % EIGHT2D
                Y_ideal = SHAPE_SIZE * sin(t);
                X_ideal = SHAPE_SIZE * sin(t) .* cos(t);
                Z_ideal = zeros(size(t));
            elseif shapeType == 4 % EIGHT3D
                S = SHAPE_SIZE;
                zScale = 0.4;
                zPhaseRad = pi/4;
                X_ideal = S .* sin(t) .* cos(t);
                Y_ideal = S .* sin(t);
                Z_ideal = (zScale * S) .* sin(t + zPhaseRad);
            else
                % Circle fallback
                X_ideal = SHAPE_SIZE * cos(t);
                Y_ideal = SHAPE_SIZE * sin(t);
                Z_ideal = zeros(size(t));
            end
            
            
            points = [X_ideal; Y_ideal; Z_ideal];
            
            if appendFlag
                % Rotate 45 deg Y (approx 0.78 rad)
                eul = [0, deg2rad(45), 0];
                R = eul2rotm(eul, 'ZYX');
                points_rot = R * points;
                points = [points, points_rot];
            end
        end
        
        function sparse_traj = sample_trajectory(dense_traj, ratio)
            sparse_traj = dense_traj(:, 1:ratio:end);
        end
    end
end