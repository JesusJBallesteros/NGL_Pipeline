% This is an example demonstrating how to correctly extract data from block
% file format files.

clc; clear;

%% ===============Set up for file parse============= %%

folder = 'C:\';
file = 'NEUR0003.DF1';
fileName = fullfile(folder, file);
fid = fopen(fileName, 'r');
rawData = uint8(fread(fid, Inf, 'uint8'));
fclose(fid);


%% ================Extract metadata from block header =============== %%


constId = (hex2num(HeaderConstants.HexConstId));
constIdBytes = typecast(constId, 'uint8');
blockStartIndices = FindDataBlockStart(rawData, constIdBytes);
numberOfBlocks = length(blockStartIndices);
startOfFirstHeader = blockStartIndices(1);
endOfFirstHeader = startOfFirstHeader + HeaderConstants.HeaderTotalBytes;
firstHeader = rawData(startOfFirstHeader:endOfFirstHeader);
HeaderStruct = ExtractHeaderData(firstHeader);

%% ================Extract neural data from blocks =============== %%

%check where each type data is in partition info
neuralIndex = find(cell2mat(arrayfun(@(x) x.DataType == uint32(DataTypeEnum.NeuralData), HeaderStruct.PartitionInfo, 'un', 0)));
isNeuralPresent = ~isempty(neuralIndex);
if isNeuralPresent
    neuralDataAsBytes = ExtractDataByType(rawData, HeaderStruct, neuralIndex, blockStartIndices, numberOfBlocks);
end


%% ================Convert neural data to physical units (V) =============== %%
% get meta data from file start event using event file reader
if isNeuralPresent
    numberOfAdcBits = 16;
    offset = 2 ^ (numberOfAdcBits - 1);
    voltageResolution = 1.95e-7;
    numberOfChannels = 64;
    frequency = 32000;

    %cast bytes to unsigned 16 bit integers and store as float
    neuralData= single(typecast(neuralDataAsBytes, 'uint16'));

    %convert uint16 values to voltage values
    for dataPointIndex = 1:length(neuralData)
        neuralData(dataPointIndex) = voltageResolution * (neuralData(dataPointIndex) - offset); 
    end

    % sort by channels
    neuralDataMat = reshape(neuralData, numberOfChannels, []);

    % get timestamps of neural data
    timestampsNeural = GetTimestamps(HeaderStruct.Timestamp, frequency, size(neuralDataMat, 2));
end
%% ================Extract audio data from blocks =============== %%

%check where each type data is in partition info
audioIndex = find(cell2mat(arrayfun(@(x) x.DataType == uint32(DataTypeEnum.Audio), HeaderStruct.PartitionInfo, 'un', 0)));
isAudioPresent = ~isempty(audioIndex);
if isAudioPresent
    audioDataAsBytes = ExtractDataByType(rawData, HeaderStruct, audioIndex, blockStartIndices, numberOfBlocks);
end

%% ================Scale audio data and save as wav =============== %%
% get meta data from file start event using event file reader

if isAudioPresent
    numberOfAudioBits = 15;
    frequency = 200000;
    isAudioSigned = true;

    audioData = ScaleAudioData(isAudioSigned, audioDataAsBytes, numberOfAudioBits); 
    % save audio data as wav file
    audiowrite('C:\Users\myPath\example200kHz.wav', audioData, frequency, 'BitsPerSample', 16);
    % get timestamps of audio data
    timestampsAudio = GetTimestamps(HeaderStruct.Timestamp, frequency, length(audioData));
end

%% ================Extract motion sensor data from blocks =============== %%

%check where each type data is in partition info
motionSensorIndex = find(cell2mat(arrayfun(@(x) x.DataType == uint32(DataTypeEnum.MotionSensor), HeaderStruct.PartitionInfo, 'un', 0)));
isMotionSensorPresent = ~isempty(motionSensorIndex);
if isMotionSensorPresent
    motionSensorAsBytes = ExtractDataByType(rawData, HeaderStruct, motionSensorIndex, blockStartIndices, numberOfBlocks);

end

%% ================Sort motion sensor data by data type =============== %%
% The values for acclMax and gyroMax are chosen by the user. They can be found using the Event
% File Viewer in the file started event.
acclMax = 2*MotionSensorConstants.G; % m/s^2, maximum value of selected range
gyroMax = 250; %  degrees per second,  maximum value of selected range
magMax = MotionSensorConstants.Magnetometer9250Range; %  Teslas,  maximum value of selected range

% extract motion sensor data from inner block

if isMotionSensorPresent

    
    motionSensorData = single(typecast(motionSensorAsBytes, 'int16'));    
    blockStartIndicesMs = FindMotionSensorBlockStart(motionSensorData, MotionSensorConstants.ConstId);    
    accelerometerDataTemp = ExtractMotionSensorDataByType(motionSensorData, blockStartIndicesMs, MotionSensorEnum.Accelerometer);
    gyroscopeDataTemp = ExtractMotionSensorDataByType(motionSensorData, blockStartIndicesMs, MotionSensorEnum.Gyroscope);
    magnetometerDataTemp = ExtractMotionSensorDataByType(motionSensorData, blockStartIndicesMs, MotionSensorEnum.Magnetometer);
    
    scaledAccelerometer = ScaleMotionSensorData(accelerometerDataTemp, MotionSensorConstants.AccelerometerNumberOfBits, acclMax);
    scaledGyroscope = ScaleMotionSensorData(gyroscopeDataTemp, MotionSensorConstants.GyroscopeNumberOfBits, gyroMax);
    scaledMagnetometer = ScaleMotionSensorData(magnetometerDataTemp, MotionSensorConstants.Magnetometer9250NumberOfBits, magMax);

    AccelerometerData = SortDataByAxis(scaledAccelerometer);
    GyroscopeData = SortDataByAxis(scaledGyroscope);
    MagnetometerData = SortDataByAxis(scaledMagnetometer);
    
    % get timestamp of a particular point
    index = 1000; % want the timestamp of the 200000th point
    timestampsAccelerometer = GetMotionSensorTimestamp(motionSensorData, blockStartIndicesMs, MotionSensorEnum.Accelerometer, MotionSensorConstants.AccelerometerFrequency);
    timestampsGyroscope = GetMotionSensorTimestamp(motionSensorData, blockStartIndicesMs, MotionSensorEnum.Gyroscope, MotionSensorConstants.GyroscopeFrequency);
    timestampsMagnetometer = GetMotionSensorTimestamp(motionSensorData, blockStartIndicesMs, MotionSensorEnum.Magnetometer, MotionSensorConstants.MagnetometerFrequency);

    % example: plot accelerometer
    plot(timestampsAccelerometer, AccelerometerData.X);
    hold on;
    plot(timestampsAccelerometer, AccelerometerData.Y);
    hold on;
    plot(timestampsAccelerometer, AccelerometerData.Z);
    ylim([-acclMax, acclMax]);
end

%% ================Extract multiple magnetometer data from blocks =============== %%

numberOfSensors = 3;
numRecordsPosition = 3; % the 4th header value is the number of records
headerLength = 4; % in 32 bit ints
recordSize = int32(36/4); %36 bytes, 4 byte ints
maxValue = 1200; %uT

%check where each type data is in partition info
multipleMagnetometerIndex = find(cell2mat(arrayfun(@(x) x.DataType == uint32(DataTypeEnum.MultipleMagnetometer), HeaderStruct.PartitionInfo, 'un', 0)));
isMultipleMagnetometerPresent = ~isempty(multipleMagnetometerIndex);  
if isMultipleMagnetometerPresent
    multipleMagnetometerAsBytes = ExtractDataByType(rawData, HeaderStruct, multipleMagnetometerIndex, blockStartIndices, numberOfBlocks);
end
blockConstId = [88 88 88 88];
blockStarts = strfind(multipleMagnetometerAsBytes.', blockConstId);
magnetometerAsInts = typecast(multipleMagnetometerAsBytes, 'int32');
blockStartInts = int32((blockStarts - 1)/4 + 1); %index of block starts for 32 bit integer data, matlab is 1 indexed

numberOfRecordsByBlock = magnetometerAsInts(blockStartInts + numRecordsPosition);

for i = 1:numberOfSensors
    Data{i}.x = zeros(sum(numberOfRecordsByBlock), 1);
    Data{i}.y = zeros(sum(numberOfRecordsByBlock), 1);
    Data{i}.z = zeros(sum(numberOfRecordsByBlock), 1);
end

startIndex = 1;
for i = 1:length(blockStartInts)
    for s = 1:numberOfSensors
        sourceStartXIdx = blockStartInts(i) + headerLength + (s-1)*numberOfSensors;
        numRecords = numberOfRecordsByBlock(i);
        Data{s}.x(startIndex:startIndex+numRecords-1) = magnetometerAsInts(sourceStartXIdx:recordSize:sourceStartXIdx + numRecords*recordSize - 1);
        Data{s}.y(startIndex:startIndex+numRecords-1) = magnetometerAsInts(sourceStartXIdx + 1:recordSize:sourceStartXIdx + numRecords*recordSize - 1);
        Data{s}.z(startIndex:startIndex+numRecords-1) = magnetometerAsInts(sourceStartXIdx + 2:recordSize:sourceStartXIdx + numRecords*recordSize - 1);
    end
            startIndex = startIndex + numRecords;

end

for s = 1:numberOfSensors
    subplot(1, numberOfSensors, s);
    plot(Data{s}.x.*10e-3) ;
    ylim([-1200 1200]);
    ylabel('\muT');
    hold on;
    plot(Data{s}.y.*10e-3);
     ylim([-1200 1200]);
    ylabel('\muT');
    hold on;
    plot(Data{s}.z.*10e-3);
     ylim([-1200 1200]);
    ylabel('\muT');
    legend('x', 'y', 'z');
end



