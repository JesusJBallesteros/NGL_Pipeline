function data = Deuteron_extractData(fid, opt)
% Deuteron_extractData  Read and parse a single Deuteron DF1 block-format file.
%
% PURPOSE:
%   Low-level reader for Deuteron DF1 files. Reads the entire file as raw
%   bytes, locates block boundaries using Deuteron block header constants,
%   extracts the requested data stream (neural / motion sensor / audio) from
%   each block, and returns scaled physical-unit data.
%   Dispatches on opt.stream to select which data type to extract.
%
% USAGE:
%   data = Deuteron_extractData(fid, opt)
%   Called from Deuteron2Kilosort and Deuteron2Fieldtrip (stream=1) and
%   from getfrom_Deuteron inside GetMotionSensors (stream=2).
%
% INPUTS:
%   fid  - file identifier from fopen on a NEUR*.DF1 file (must be open)
%   opt  - struct with fields:
%            .stream     integer: 1=neural, 2=motion sensor, 3=audio
%            .acclMax    accelerometer full-scale range (m/s²)   [stream=2 only]
%            .gyroMax    gyroscope full-scale range (deg/s)      [stream=2 only]
%            .magMax     magnetometer full-scale range (Tesla)   [stream=2 only]
%
% OUTPUTS:
%   stream=1 (neural):
%     data   - [1 × nSamples single] raw ADC values (uint16 cast to single).
%              Caller applies ADC→µV conversion after this call.
%   stream=2 (motion sensor):
%     data   - struct with fields:
%                .Accelerometer.Data  [N × 3 single] in m/s²
%                .Gyroscope.Data      [N × 3 single] in deg/s
%                .Magnetometer.Data   [N × 3 single] in Tesla
%                .<sensor>.timestamps [N × 1 double] in ms
%
% CALLS:
%   Deuteron toolbox: HeaderConstants, FindDataBlockStart, ExtractHeaderData,
%   ExtractDataByType, FindMotionSensorBlockStart, ExtractMotionSensorDataByType,
%   ScaleMotionSensorData, SortDataByAxis, GetMotionSensorTimestamp,
%   MotionSensorConstants, MotionSensorEnum, DataTypeEnum
%
% Last modified 08.05.2026 (Jesus)

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