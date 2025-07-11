function [time, nSamples] = intan_gettimeseries(sessions, ss)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

tmp.fileinfo = dir('time.dat'); % Name won't change
tmp.num_samples = tmp.fileinfo.bytes/4; % int32 = 4 bytes
tmp.samplerate = sessions.info{1,ss}.amplifier_sample_rate;

% Open file according to INTAN
fid = fopen(tmp.fileinfo.name, 'r');

% Extract time data (sample times). It comes in samples, we 
% convert to seconds based on sample rate. 
time = fread(fid, tmp.num_samples, 'int32') / tmp.samplerate;
nSamples = tmp.num_samples;

% Close file.
fclose(fid); clear tmp fid ans

end