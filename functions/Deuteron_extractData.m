function data = Deuteron_extractData(fid, opt)
% This extracts data from Deuteron's block file format files.
% Determines how to proceed depending on the value of the 'stream' variable.

% Jesus 05.01.2024. Modified from Deuteron's

%% Parse file
rawData = uint8(fread(fid, Inf, 'uint8'));

%% Extract metadata from block header 
constId             = (hex2num(HeaderConstants.HexConstId));
constIdBytes        = typecast(constId, 'uint8');
blockStartIndices   = FindDataBlockStart(rawData, constIdBytes);
numberOfBlocks      = length(blockStartIndices);
startOfFirstHeader  = blockStartIndices(1);
endOfFirstHeader    = startOfFirstHeader + HeaderConstants.HeaderTotalBytes;
firstHeader         = rawData(startOfFirstHeader:endOfFirstHeader);
HeaderStruct        = ExtractHeaderData(firstHeader);

switch opt.stream
    case 1
    %% Extract neural data from blocks. 
    % Check where each type data is in partition info.
    neuralIndex = find(cell2mat(arrayfun(@(x) x.DataType == uint32(DataTypeEnum.NeuralData), HeaderStruct.PartitionInfo, 'un', 0)));
    isNeuralPresent = ~isempty(neuralIndex);
    if isNeuralPresent
        neuralDataAsBytes = ExtractDataByType(rawData, HeaderStruct, neuralIndex, blockStartIndices, numberOfBlocks);
        
        % Cast bytes to unsigned 16 bit integers and store as float
        data = single(typecast(neuralDataAsBytes, 'uint16'));

        % get timestamps of neural data
        %    out.timestampsNeural = GetTimestamps(HeaderStruct.Timestamp, frequency, size(neuralData, 2));
    end
    
    case 2
    %% Extract motion sensor data from blocks.
    % Check where each type data is in partition info.
    motionSensorIndex = find(cell2mat(arrayfun(@(x) x.DataType == uint32(DataTypeEnum.MotionSensor), HeaderStruct.PartitionInfo, 'un', 0)));
    isMotionSensorPresent = ~isempty(motionSensorIndex);
    if isMotionSensorPresent
        motionSensorAsBytes = ExtractDataByType(rawData, HeaderStruct, motionSensorIndex, blockStartIndices, numberOfBlocks);

        % extract motion sensor data from inner block    
        motionSensorData    = single(typecast(motionSensorAsBytes, 'int16'));    
        blockStartIndicesMs = FindMotionSensorBlockStart(motionSensorData, MotionSensorConstants.ConstId);    
        accelerometerDataTemp = ExtractMotionSensorDataByType(motionSensorData, blockStartIndicesMs, MotionSensorEnum.Accelerometer);
        gyroscopeDataTemp   = ExtractMotionSensorDataByType(motionSensorData, blockStartIndicesMs, MotionSensorEnum.Gyroscope);
        magnetometerDataTemp = ExtractMotionSensorDataByType(motionSensorData, blockStartIndicesMs, MotionSensorEnum.Magnetometer);
        
        scaledAccelerometer = ScaleMotionSensorData(accelerometerDataTemp, MotionSensorConstants.AccelerometerNumberOfBits, opt.acclMax);
        scaledGyroscope     = ScaleMotionSensorData(gyroscopeDataTemp, MotionSensorConstants.GyroscopeNumberOfBits, opt.gyroMax);
        scaledMagnetometer  = ScaleMotionSensorData(magnetometerDataTemp, MotionSensorConstants.Magnetometer9250NumberOfBits, opt.magMax);
    
        data.Accelerometer.Data  = SortDataByAxis(scaledAccelerometer);
        data.Gyroscope.Data      = SortDataByAxis(scaledGyroscope);
        data.Magnetometer.Data   = SortDataByAxis(scaledMagnetometer);
        
        % get timestamp of a particular point
        data.Accelerometer.timestamps   = GetMotionSensorTimestamp(motionSensorData, blockStartIndicesMs, MotionSensorEnum.Accelerometer, MotionSensorConstants.AccelerometerFrequency);
        data.Gyroscope.timestamps       = GetMotionSensorTimestamp(motionSensorData, blockStartIndicesMs, MotionSensorEnum.Gyroscope, MotionSensorConstants.GyroscopeFrequency);
        data.Magnetometer.timestamps    = GetMotionSensorTimestamp(motionSensorData, blockStartIndicesMs, MotionSensorEnum.Magnetometer, MotionSensorConstants.MagnetometerFrequency);
    end

    case 3 
    %% Extract audio data from blocks
    % Check where each type data is in partition info
    audioIndex = find(cell2mat(arrayfun(@(x) x.DataType == uint32(DataTypeEnum.Audio), HeaderStruct.PartitionInfo, 'un', 0)));
    isAudioPresent = ~isempty(audioIndex);
    if isAudioPresent
        audioDataAsBytes = ExtractDataByType(rawData, HeaderStruct, audioIndex, blockStartIndices, numberOfBlocks);
    
        % Scale audio data and save as wav
        % Get meta data from file start event using event file reader
        numberOfAudioBits = 15;
        frequency = 200000;
        isAudioSigned = true;
    
        data.audioData = ScaleAudioData(isAudioSigned, audioDataAsBytes, numberOfAudioBits); 
        
        % get timestamps of audio data
        data.timestampsAudio = GetTimestamps(HeaderStruct.Timestamp, frequency, length(audioData));

        % % save audio data as wav file bring out
        % audiowrite('C:\Users\myPath\example200kHz.wav', audioData, frequency, 'BitsPerSample', 16);
     end

    case 4
    %% Reserved
end

end