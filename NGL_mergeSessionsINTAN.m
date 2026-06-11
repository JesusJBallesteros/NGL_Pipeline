%% NGL_mergeSessionsINTAN. Concatenate two INTAN sessions into one.
%
% PURPOSE:
%   One-off script for the special case where two INTAN recordings of the
%   SAME subject (and same probe / sample rate / channel count) should be
%   treated downstream as a single continuous session — e.g. when a single
%   experimental day was split into two recording files by accident, or
%   when a "before / after" manipulation produced two halves that need to
%   be sorted together so cluster identity stays consistent across them.
%
%   Produces most of the per-session artefacts NGL01_Main would have made,
%   EXCEPT it does NOT run Kilosort, Bombcell or write NWB. Run NGL01_Main
%   (with opt.kilosort=true, opt.bombcell=true) against the merged folder
%   afterwards to sort and curate; clusters then span the whole [A B]
%   recording with a single identity.
%
% USAGE (from NGL_SetAndRunMe):
%   subject     = 'JN3';                          % subject ID (string or char)
%   mergeDates  = {'20260101' '20260103'};        % {A, B} in chronological order
%   mergeTag    = 'JN3_pre_post';                 % user-supplied output folder name
%   % standard inputs as in any NGL run:
%   datadrive   = 'D';
%   studyname   = 'colorWheel';
%   opt         = struct();
%   opt.numChannels   = 32;
%   opt.RetrieveEvents = true;
%   opt.alignto       = {'itiOn'};
%   opt.GetMotionSensors = true;                  % if AUX has motion data
%   opt.FieldTrip     = true;
%   % opt.kilosort / .bombcell are forced false by this script.
%   NGL_mergeSessionsINTAN
%
% INPUTS (workspace variables):
%   subject     - char/string, exactly ONE subject ID.
%   mergeDates  - 1x2 cell of 'YYYYMMDD' char vectors. The first becomes
%                 session A (samples 1..nA in the merged record); the
%                 second becomes session B (samples nA+1..nA+nB).
%   mergeTag    - char, output folder name under
%                 <datadrive>:\<studyname>\data\preprocessing\<subject>\.
%                 Pick a name that's unmistakeably NOT a date (e.g.
%                 'JN3_pre_post' or 'D14_AB').
%   opt         - standard NGL options struct. Used fields: numChannels,
%                 alignto, RetrieveEvents, GetMotionSensors, FieldTrip,
%                 lowpass, lowpassFT, highpass, linefilter, CAR,
%                 dwnsmplRate, lowpass, lowpassFT, addtime, uselog.
%                 The script FORCES opt.kilosort = opt.bombcell = false
%                 and opt.doNWB = false so the heavy stages are skipped.
%
% OUTPUT (in <input.processed>/<subject>/<mergeTag>/):
%   <mergeTag>.bin           concatenated raw [A B] int16 [nCh x (nA+nB)]
%   EventRecord.mat          merged digital events; B timestamps shifted
%                            by nSamplesA; TimeBreak marks the boundary.
%   MotionData_raw.mat       merged AUX motion samples (if GetMotionSensors).
%   <mergeTag>_FTcont.mat    merged FieldTrip continuous data.
%   events.mat / trialdef.mat re-derived from merged EventRecord via
%                            trialdefGen so downstream NGL02 paths work.
%   mergeMeta.mat            metadata for downstream: sessionA, sessionB,
%                            sampleBoundary, nTrialsA, nTrialsB,
%                            trialSessionTag, mergeTimestamp.
%   mergeManifest.txt        human-readable summary of the merge.
%   preprocInfo.mat          snapshot of the opt that was used.
%   staging/<dateA>/, staging/<dateB>/   per-session intermediates
%                            (kept for inspection; delete manually or set
%                             opt.merge.cleanupStaging = true to remove).
%
% PIPELINE:
%   00. Prep + validate inputs
%   01. Per-session extraction into staging/{A,B}/
%   02. Compute sample boundary nSamplesA from staged .bin file size
%   03. Concatenate analog (.bin), digital (EventRecord), AUX (MotionData_raw),
%       and FieldTrip continuous; mark TimeBreak at the join.
%   04. Re-derive events.mat + trialdef.mat from merged EventRecord.
%   05. Write mergeMeta + preprocInfo + manifest.
%   06. Optional staging cleanup.
%
% RUN KILOSORT / BOMBCELL AFTERWARDS:
%   In NGL_SetAndRunMe:
%     subjects = {'JN3'};
%     dates    = {'JN3_pre_post'};   % the merge tag, treated as a session
%     opt.kilosort = true; opt.bombcell = true;
%     NGL01_Main
%
% LIMITATIONS:
%   - INTAN only (asserts both sessions report sessions(x).info.fileformat
%     in {'fileperch','filepertype','tradFormat'}).
%   - Channel count and sample rate must match between A and B.
%   - NWB export is not run; regenerate from the merged .bin if needed.
%   - Wall-clock timestamps (EventRecord.TimeMsFromMidnight) are preserved
%     verbatim per session (NOT shifted) since they reflect real-world
%     clock; only sample-indexed timestamps are shifted.
%
% Last modified 09.06.2026 (Jesus)

%% 00. Standard scaffolding + input validation.
NGL00_Prep

assert(exist('subject','var') == 1 && (ischar(subject) || isstring(subject)) && ~isempty(subject), ...
    'NGL:mergeINTAN', 'Workspace variable `subject` (single subject ID) is required.');
subject = char(subject);

assert(exist('mergeDates','var') == 1 && iscell(mergeDates) && numel(mergeDates) == 2 && ...
       all(cellfun(@(d) ischar(d) && length(d) == 8, mergeDates)), ...
    'NGL:mergeINTAN', ...
    'Workspace variable `mergeDates` must be a 1x2 cell of ''YYYYMMDD'' chars.');
dateA = mergeDates{1};
dateB = mergeDates{2};

assert(exist('mergeTag','var') == 1 && ischar(mergeTag) && ~isempty(mergeTag), ...
    'NGL:mergeINTAN', ...
    'Workspace variable `mergeTag` must be a non-empty char (the merged session''s folder name).');
assert(length(mergeTag) ~= 8 || any(~isstrprop(mergeTag, 'digit')), ...
    'NGL:mergeINTAN', ...
    ['mergeTag ''%s'' looks like an 8-digit date. Choose a non-date name ', ...
     '(e.g. ''%s_pre_post'') so downstream code can''t mistake it for a real recording.'], ...
    mergeTag, subject);

% Force heavy stages OFF for the merger.
opt.kilosort  = false;
opt.bombcell  = false;
opt.callBcGUI = false;
opt.phy       = false;
opt.doNWB     = false;

% Constrain the discovery scope: this subject, these two dates.
input.subjects = {subject};
input.dates    = mergeDates;

[input, opt] = set_default(input, opt);
input.sessions = findSessions(input);

assert(input.nsubjects == 1, ...
    'NGL:mergeINTAN', 'Expected exactly 1 subject; got %d.', input.nsubjects);
assert(input.sessions(1).nsessions == 2, ...
    'NGL:mergeINTAN', ...
    'Expected 2 sessions matching mergeDates; found %d under %s.', ...
    input.sessions(1).nsessions, input.sessions(1).folder);

% Confirm the two found sessions match mergeDates in order.
foundDates = input.sessions(1).list;
assert(any(strcmp(foundDates, dateA)) && any(strcmp(foundDates, dateB)), ...
    'NGL:mergeINTAN', ...
    'mergeDates {%s, %s} not found in discovered sessions {%s}.', ...
    dateA, dateB, strjoin(foundDates, ', '));
% Force ordering A then B regardless of what findSessions returned.
idxA = find(strcmp(foundDates, dateA), 1);
idxB = find(strcmp(foundDates, dateB), 1);
input.sessions(1).list = foundDates([idxA, idxB]);

%% 01. Per-session extraction into staging/{A,B}/.
mergedFolder = fullfile(input.processed, subject, mergeTag);
stagingDir   = fullfile(mergedFolder, 'staging');
if ~exist(mergedFolder, 'dir'), mkdir(mergedFolder); end
if ~exist(stagingDir,   'dir'), mkdir(stagingDir);   end

stagedInfo = struct('date', {dateA, dateB}, ...
                    'stage', {fullfile(stagingDir, dateA), fullfile(stagingDir, dateB)}, ...
                    'binFile', {'', ''}, ...
                    'eventFile', {'', ''}, ...
                    'motionFile', {'', ''}, ...
                    'ftFile', {'', ''}, ...
                    'nChannels', {NaN, NaN}, ...
                    'sampleRate', {NaN, NaN}, ...
                    'nSamples', {NaN, NaN}, ...
                    'nTrials', {NaN, NaN});

for r = 1:2
    if ~exist(stagedInfo(r).stage, 'dir'), mkdir(stagedInfo(r).stage); end
    input.run = [1, r];     % (subject index, session index) per NGL convention

    % Per-session opt scaffolding. prepforsession sets opt.FolderProcDataMat
    % to <input.processed>/<subject>/<date>/; we'll move outputs out of that
    % folder into our staging dir at the end of this iteration to avoid
    % polluting downstream NGL01 runs.
    [input, opt] = prepforsession(input, opt);
    sessOutDir   = opt.FolderProcDataMat;
    if ~exist(sessOutDir, 'dir'), mkdir(sessOutDir); end

    sessMeta = input.sessions(1).info;
    assert(ismember(sessMeta.fileformat, {'fileperch','filepertype','tradFormat'}), ...
        'NGL:mergeINTAN', ...
        'Session %s reports fileformat ''%s'' which is not an INTAN format. The merger is INTAN-only.', ...
        stagedInfo(r).date, sessMeta.fileformat);

    stagedInfo(r).nChannels  = sessMeta.nChannels;
    stagedInfo(r).sampleRate = sessMeta.amplifier_sample_rate;

    fprintf('NGL_mergeSessionsINTAN: extracting session %s (%s, %d ch @ %g Hz)\n', ...
            stagedInfo(r).date, sessMeta.fileformat, sessMeta.nChannels, ...
            sessMeta.amplifier_sample_rate);

    % 01.a  .bin file (analog Kilosort-ready).
    binDst = fullfile(stagedInfo(r).stage, [stagedInfo(r).date '.bin']);
    if ~isfile(binDst)
        Intan2Kilosort_wrapper(input.sessions(1), opt);
        binSrc = fullfile(sessOutDir, [stagedInfo(r).date '.bin']);
        assert(isfile(binSrc), 'NGL:mergeINTAN', ...
            'Intan2Kilosort_wrapper did not produce %s', binSrc);
        binDst = fullfile(stagedInfo(r).stage, [stagedInfo(r).date '.bin']);
        movefile(binSrc, binDst);
        stagedInfo(r).binFile  = binDst;
        fInfo                  = dir(binDst);
        stagedInfo(r).nSamples = fInfo.bytes / (2 * sessMeta.nChannels);
    else
        fprintf('NGL_mergeSessionsINTAN: a file %s already exists.\n', binDst);
        stagedInfo(r).binFile  = binDst;
        fInfo                  = dir(binDst);
        stagedInfo(r).nSamples = fInfo.bytes / (2 * sessMeta.nChannels);
    end

    % 01.b  Events.
    evtDst = fullfile(stagedInfo(r).stage, 'EventRecord.mat');
    if opt.RetrieveEvents && ~isfile(evtDst)
        EventRecord = INTAN_ExtractEvents(input, opt);
        save(evtDst, 'EventRecord');
        stagedInfo(r).eventFile = evtDst;

        % Read eventcode list
        if ~isfield(opt,'eventdef')       
            opt.eventdef  = eventDefinitions(input.sessions(input.run(1)).info.fileformat);
        end 

        % Count trials in this session for the merge metadata.
        try
            tmpEvents = events2align(opt);
            [~, tmpTrialdef, ~, ~] = trialdefGen(EventRecord, opt);
            stagedInfo(r).nTrials = size(tmpTrialdef{2,1}, 1);
        catch ME
            warning('NGL:mergeINTAN:trialCount', ...
                'Could not count trials for %s (%s); leaving nTrials=NaN. Will renumber post-merge anyway.', ...
                stagedInfo(r).date, ME.message);
        end
    end

    % 01.c  AUX motion data.
    motDst = fullfile(stagedInfo(r).stage, 'MotionData_raw.mat');
    if opt.GetMotionSensors && ~isfile(motDst)
        % INTAN motion sensors are extracted into MotionData_raw.mat
        try
            GetMotionSensors(opt, input);
            motSrc = fullfile(sessOutDir, 'MotionData_raw.mat');
            if isfile(motSrc)
                motDst = fullfile(stagedInfo(r).stage, 'MotionData_raw.mat');
                movefile(motSrc, motDst);
                stagedInfo(r).motionFile = motDst;
            end
        catch ME
            warning('NGL:mergeINTAN:motion', ...
                'ProcessMotionSensors failed for %s (%s); merged record will lack MotionData_raw.', ...
                stagedInfo(r).date, ME.message);
        end
    end

    % 01.d  FieldTrip continuous.
    ftDst = fullfile(stagedInfo(r).stage, [stagedInfo(r).date '_FTcont.mat']);
    if opt.FieldTrip && ~isfile(ftDst)
        try
            INTANdata = intan2MAT_wrapper(input.sessions(input.run(1)), opt);
            MAT2FieldTrip(INTANdata, opt);
            ftSrc = fullfile(sessOutDir, [stagedInfo(r).date '_FTcont.mat']);
            if ~isfile(ftSrc)
                warning('NGL:mergeINTAN:noFT', ...
                    ['FT continuous file %s not found. Re-run NGL01_Main per session ', ...
                     'to produce *_FTcont.mat, then re-run the merger.'], ftSrc);
            else
                ftDst = fullfile(stagedInfo(r).stage, [stagedInfo(r).date '_FTcont.mat']);
                movefile(ftSrc, ftDst);
                stagedInfo(r).ftFile = ftDst;
            end
        catch ME
            warning('NGL:mergeINTAN:FT', 'FT continuous extraction failed for %s (%s).', ...
                stagedInfo(r).date, ME.message);
        end
    elseif isfile(ftDst)
        fprintf('NGL_mergeSessionsINTAN: a FTfile %s already exists.\n', ftDst);
    end

    % Tidy up: drop the now-empty per-session NGL01 output dir.
    if exist(sessOutDir, 'dir') && numel(dir(fullfile(sessOutDir, '*')))==2
        rmdir(sessOutDir, 's');
    end
end

%% 02. Boundary checks across A and B.
assert(stagedInfo(1).nChannels == stagedInfo(2).nChannels, ...
    'NGL:mergeINTAN', ...
    'Channel count mismatch: %d (%s) vs %d (%s). Cannot merge.', ...
    stagedInfo(1).nChannels, dateA, stagedInfo(2).nChannels, dateB);
assert(stagedInfo(1).sampleRate == stagedInfo(2).sampleRate, ...
    'NGL:mergeINTAN', ...
    'Sample-rate mismatch: %g Hz (%s) vs %g Hz (%s). Cannot merge.', ...
    stagedInfo(1).sampleRate, dateA, stagedInfo(2).sampleRate, dateB);

nChannels  = stagedInfo(1).nChannels;
sampleRate = stagedInfo(1).sampleRate;
nSamplesA  = stagedInfo(1).nSamples;
nSamplesB  = stagedInfo(2).nSamples;

fprintf('NGL_mergeSessionsINTAN: boundary at sample %d (%g s); merged total %d samples.\n', ...
        nSamplesA, nSamplesA / sampleRate, nSamplesA + nSamplesB);

%% 03. Concatenate.
% 03.a  .bin file: stream-copy A then B into <mergeTag>.bin.
mergedBin = fullfile(mergedFolder, [mergeTag '.bin']);
fprintf('NGL_mergeSessionsINTAN: writing %s\n', mergedBin);
if ~isfile(mergedBin)
    concatBinFiles(stagedInfo(1).binFile, stagedInfo(2).binFile, mergedBin);
else
    fprintf('NGL_mergeSessionsINTAN: concatenated .bin file %s already exists.\n', mergedBin);
end

% 03.b  EventRecord: shift B sample timestamps by nSamplesA, renumber,
%       set TimeBreak to the boundary.
%       Re-derive events.mat and trialdef.mat from the merged EventRecord.
if isfile(fullfile(mergedFolder, 'EventRecord.mat'))
    disp('NGL_mergeSessionsINTAN: concatenated Event files already exists.\n');
else
    if ~isempty(stagedInfo(1).eventFile) && ~isempty(stagedInfo(2).eventFile)
        erA = load(stagedInfo(1).eventFile, 'EventRecord');
        erB = load(stagedInfo(2).eventFile, 'EventRecord');
        erA = erA.EventRecord;
        erB = erB.EventRecord;
    
        EventRecord = struct();
        EventRecord.EventType   = [erA.EventType(:);   erB.EventType(:)];
        EventRecord.TimeStamp   = [erA.TimeStamp(:);   erB.TimeStamp(:) + nSamplesA];
        EventRecord.EventNumber = (1:numel(EventRecord.EventType))';
        if isfield(erA, 'TimeMsFromMidnight')
            % Wall-clock times preserved verbatim per session; downstream that
            % uses these must be aware of the discontinuity.
            EventRecord.TimeMsFromMidnight = [erA.TimeMsFromMidnight(:); erB.TimeMsFromMidnight(:)];
        end
        if isfield(erA, 'TimeSource')
            EventRecord.TimeSource = [erA.TimeSource(:); erB.TimeSource(:)];
        end
        if isfield(erA, 'Details')
            EventRecord.Details    = [erA.Details(:);    erB.Details(:)];
        end
        % Mark the A->B join as the (only) time break in the merged record.
        EventRecord.TimeBreak  = {numel(erA.EventType),numel(erA.EventType)+1};
        clear erA erB
        
        % Re-derive events.mat and trialdef.mat from the merged EventRecord.
        [events, trialdef, eventdef, EventRecord] = trialdefGen(EventRecord, opt);
        save(fullfile(mergedFolder, 'EventRecord.mat'), 'EventRecord');
        save(fullfile(mergedFolder, 'events.mat'), 'events');
        save(fullfile(mergedFolder, 'trialdef.mat'), 'trialdef');
    end
end

% 03.c  MotionData_raw: concat per-channel samples.
if isfile(fullfile(mergedFolder, 'MotionData_raw.mat'))
    disp('NGL_mergeSessionsINTAN: MotionData concatenated file already exists.\n');
else
    if ~isempty(stagedInfo(1).motionFile) && ~isempty(stagedInfo(2).motionFile)
        motA = load(stagedInfo(1).motionFile);
        motB = load(stagedInfo(2).motionFile);
        motA = motA.raw;
        motB = motB.raw;
        MotionData_raw = concatMotionRaw(motA, motB);
        save(fullfile(mergedFolder, 'MotionData_raw.mat'), 'MotionData_raw');
        clear motA motB MotionData_raw
    end
end

% 03.d  FieldTrip continuous: horzcat trial{1}, shift time{1}.
if isfile(fullfile(mergedFolder, [mergeTag '_FTcont.mat']))
    disp('NGL_mergeSessionsINTAN: FT concatenated file already exists.\n');
else
    if ~isempty(stagedInfo(1).ftFile) && ~isempty(stagedInfo(2).ftFile)
        SA = load(stagedInfo(1).ftFile);
        SB = load(stagedInfo(2).ftFile);
        % FT files usually save FT_data; tolerate either name.
        if isfield(SA, 'FT_data'), ftA = SA.FT_data; else, ftA = struct2first(SA); end
        if isfield(SB, 'FT_data'), ftB = SB.FT_data; else, ftB = struct2first(SB); end
        FT_data = concatFTcontinuous(ftA, ftB);
        save(fullfile(mergedFolder, [mergeTag '_FTcont.mat']), 'FT_data', '-v7.3');
        clear SA SB ftA ftB FT_data
    end
end

%% 04. Write mergeMeta + preprocInfo + manifest.
nTrialsA = stagedInfo(1).nTrials;
nTrialsB = stagedInfo(2).nTrials;
if isnan(nTrialsA), nTrialsA = 0; end
if isnan(nTrialsB), nTrialsB = 0; end

mergeMeta = struct();
mergeMeta.sessionA          = dateA;
mergeMeta.sessionB          = dateB;
mergeMeta.subject           = subject;
mergeMeta.mergeTag          = mergeTag;
mergeMeta.sampleBoundary    = nSamplesA;
mergeMeta.sampleRate        = sampleRate;
mergeMeta.nChannels         = nChannels;
mergeMeta.nSamplesA         = nSamplesA;
mergeMeta.nSamplesB         = nSamplesB;
mergeMeta.nSamplesMerged    = nSamplesA + nSamplesB;
mergeMeta.nTrialsA          = nTrialsA;
mergeMeta.nTrialsB          = nTrialsB;
mergeMeta.trialBoundary     = nTrialsA;    % trial index at which B starts (i.e. trial nA+1 is first of B)
mergeMeta.trialSessionTag   = [repmat({'A'}, nTrialsA, 1); repmat({'B'}, nTrialsB, 1)];
mergeMeta.mergeTimestamp    = datetime('now');
mergeMeta.mergerVersion     = '1.0';
save(fullfile(mergedFolder, 'mergeMeta.mat'), 'mergeMeta');

% preprocInfo snapshot (re-use the same helper NGL01 uses).
try
    savePreprocInfo(opt, input, mergedFolder);
catch ME
    warning('NGL:mergeINTAN:preprocInfo', ...
        'Could not save preprocInfo snapshot (%s); writing a minimal fallback.', ME.message);
    preprocInfo = struct('opt', opt, 'input', input, 'mergeMeta', mergeMeta); %#ok<NASGU>
    save(fullfile(mergedFolder, 'preprocInfo.mat'), 'preprocInfo');
end

% Manifest (human-readable).
manifestPath = fullfile(mergedFolder, 'mergeManifest.txt');
fid = fopen(manifestPath, 'w');
fprintf(fid, 'NGL_mergeSessionsINTAN  manifest\n');
fprintf(fid, '================================\n');
fprintf(fid, 'Merge timestamp:  %s\n', char(mergeMeta.mergeTimestamp));
fprintf(fid, 'Subject:          %s\n', subject);
fprintf(fid, 'Merge tag:        %s\n', mergeTag);
fprintf(fid, 'Session A:        %s  (%d samples, %d trials)\n', dateA, nSamplesA, nTrialsA);
fprintf(fid, 'Session B:        %s  (%d samples, %d trials)\n', dateB, nSamplesB, nTrialsB);
fprintf(fid, 'Channels:         %d  @  %g Hz\n', nChannels, sampleRate);
fprintf(fid, 'Sample boundary:  %d  (= %.3f s into recording)\n', nSamplesA, nSamplesA/sampleRate);
fprintf(fid, 'Trial boundary:   trial %d (B starts at trial %d in the merged sequence)\n', ...
        nTrialsA, nTrialsA + 1);
fprintf(fid, '\nDownstream:\n');
fprintf(fid, '  - Run NGL01_Main with dates = {''%s''} to spike-sort and curate.\n', mergeTag);
fprintf(fid, '  - condition_script(_) or NGL02_postPhy should call attachMergeSession\n');
fprintf(fid, '    to add condition.session = {''A'', ''A'', ..., ''B'', ''B'', ...} per trial.\n');
fclose(fid);
fprintf('NGL_mergeSessionsINTAN: %s\n', manifestPath);

% %% 06. Optional staging cleanup.
% if isfield(opt,'merge') && isfield(opt.merge,'cleanupStaging') && opt.merge.cleanupStaging
%     rmdir(stagingDir, 's');
%     fprintf('NGL_mergeSessionsINTAN: staging cleaned.\n');
% else
%     fprintf('NGL_mergeSessionsINTAN: staging kept at %s (delete manually or set opt.merge.cleanupStaging=true).\n', ...
%             stagingDir);
% end

% --- Local UI helpers ----------------------------------------------------

function s = struct2first(S)
% Return the value of the first non-empty field of S — used as a tolerant
% fallback when we don't know whether the FT file saved its struct under
% 'FT_data' or under some other name.
    fns = fieldnames(S);
    for k = 1:numel(fns)
        if ~isempty(S.(fns{k}))
            s = S.(fns{k});
            return
        end
    end
    s = struct();
end
