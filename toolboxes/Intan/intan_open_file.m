% Open file and give a UI status warning in case of failure
function [success, fid] = intan_open_file(filename)
[fid, errmsg] = fopen(filename, 'r');
if ~strcmp(errmsg, '')
    status_str = append('Failed to read ', filename, '. Is it in this directory?');
    update_ui_status(status_str, 'red');
    success = false;
    return;
else
    success = true;
end
return;
end
