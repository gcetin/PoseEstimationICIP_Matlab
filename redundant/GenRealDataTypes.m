classdef GenRealDataTypes
    properties (Constant)
        MAX_IMAGE_WIDTH = 1920;
        MAX_IMAGE_HEIGHT = 1080;
        F_mm = 3.3333333333333333333333 * 1.5;
        Px_mm = 6.4;
        Py_mm = 3.6;
        
        fx_ = 1.000621893257644*1.0e+03;
        fy_ = 1.017473951321606*1.0e+03;
        ux_ = 1.0e+03 * 0.968620789005673;
        uy_ = 1.0e+03 * 0.538724895794938;

        MIN_Z_DIST_FROM_CAM = 15.0;
        
        % Shape Types (Enums simulated as constants)
        SHAPE_EIGHT2D = 0;
        SHAPE_QUATREFOIL = 1;
        SHAPE_CIRCLE = 2;
        SHAPE_HEART = 3;
        SHAPE_EIGHT3D = 4;
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
            val = obj.fx_;
        end
        function val = get.fy(obj)
             val = obj.fy_;
        end
        function val = get.ux(obj)
            val = obj.ux_;
        end
        function val = get.uy(obj)
            val = obj.uy_;
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