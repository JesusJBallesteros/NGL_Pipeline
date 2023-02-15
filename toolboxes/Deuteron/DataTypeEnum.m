classdef DataTypeEnum < uint32
    %DATATYPEENUM Summary of this class goes here
    %   Detailed explanation goes here
    
    enumeration
        NoData (0),
        EventData (1),
        NeuralData (2),
        MotionSensor (3),
        Audio (4),
        MultipleMagnetometer (8)
    end
        
end

