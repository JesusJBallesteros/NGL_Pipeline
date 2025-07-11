 function DisplayRecord( record )
%DISPLAYRECORDASSTRING takes one input, a record, and writes the
%information contained in the record to the command window.
    fprintf('Event number: %d\nEventType: %s\nTimeStamp: %s\nTimeMsFromMidnight: %d\nTimeSource: %s\nDetails: %s\n', ...
        str2double(char(record(1))), char(record(2)), char(record(3)), str2double(char(record(4))), char(record(5)), ...
        char(record(6))); %display chosen record
end

