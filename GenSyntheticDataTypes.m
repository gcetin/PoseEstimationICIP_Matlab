classdef GenSyntheticDataTypes
    properties (Constant)
        MAX_IMAGE_WIDTH = 1920;
        MAX_IMAGE_HEIGHT = 1080;
        F_mm = 3.3333333333333333333333 * 1.5;
        Px_mm = 6.4;
        Py_mm = 3.6;
        
        MIN_Z_DIST_FROM_CAM = 15.0;
        
        % Shape Types (Enums simulated as constants)
        SHAPE_EIGHT = 0;
        SHAPE_QUATREFOIL = 1;
        SHAPE_CIRCLE = 2;
        SHAPE_HEART = 3;
    end
    
    properties (Dependent)
        fx
        fy
        ux
        uy
        KK
        distCoeffs
        deltaT
    end
    
    methods
        function val = get.fx(obj)
            val = obj.F_mm / (obj.Px_mm / obj.MAX_IMAGE_WIDTH);
        end
        function val = get.fy(obj)
             val = obj.F_mm / (obj.Py_mm / obj.MAX_IMAGE_HEIGHT);
        end
        function val = get.ux(obj)
            val = obj.MAX_IMAGE_WIDTH / 2.0;
        end
        function val = get.uy(obj)
            val = obj.MAX_IMAGE_HEIGHT / 2.0;
        end
        function val = get.KK(obj)
            val = eye(3);
            val(1,1) = obj.fx;
            val(2,2) = obj.fy;
            val(1,3) = obj.ux;
            val(2,3) = obj.uy;
        end
        function val = get.distCoeffs(obj)
            val = zeros(1,5);
        end
        function val = get.deltaT(obj)
            val = 0.005;
        end
    end
end