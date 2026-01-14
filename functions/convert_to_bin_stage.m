function [binPath, didConvert] = convert_to_bin_stage(input, opt)
%CONVERT_TO_BIN_STAGE Suggested shared conversion stage for .bin creation.
% This helper centralizes the "high-pass to .bin" step for Deuteron and
% Intan formats. It is NOT wired into the pipeline yet; call it explicitly
% from wrappers once validated in your environment.
%
% INPUTS:
%   input : pipeline input struct (from set_default/prepforsession)
%   opt   : options struct (must include opt.ext, opt.bin, opt.FolderProcDataMat)
%
% OUTPUTS:
%   binPath    : full path to the .bin file (if present or created)
%   didConvert : true if this function created a new .bin file
%
% Example usage (Deuteron):
%   [binPath, didConvert] = convert_to_bin_stage(input, opt);
%
% Example usage (Intan):
%   opt.ext = input.sessions(input.run(1)).info.fileformat;
%   [binPath, didConvert] = convert_to_bin_stage(input, opt);
%
% Jesus / NGL pipeline - placeholder abstraction

didConvert = false;
binPath = fullfile(opt.FolderProcDataMat, [opt.SavFileName, '.bin']);

if ~isfield(opt, 'bin') || ~opt.bin
    return;
end

if isfile(binPath)
    return;
end

switch opt.ext
    case {'DT2', 'DF1'}
        Deuteron2Kilosort(opt);
        didConvert = isfile(binPath);
    case {'fileperch', 'filepertype', 'tradFormat'}
        Intan2Kilosort_wrapper(input.sessions(input.run(1)), opt);
        didConvert = isfile(binPath);
    otherwise
        warning('convert_to_bin_stage:unsupportedFormat', ...
            'Unsupported format for bin conversion: %s', opt.ext);
end
end
