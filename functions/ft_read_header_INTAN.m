function [hdr] = ft_read_header_INTAN(filename, varargin)

% ft_read_header MODIFICATION BY JESUS on 26.09.2022 from original
% reads header information from a variety of EEG, MEG and other time
% series data files and represents the header information in a common
% data-independent structure. The supported formats are listed below.
%
% Use as
%   hdr = ft_read_header(filename, ...)
%
% Additional options should be specified in key-value pairs and can be
%   'headerformat'   = string
%   'fallback'       = can be empty or 'biosig' (default = [])
%   'checkmaxfilter' = boolean, whether to check that maxfilter has been correctly applied (default = true)
%   'chanindx'       = list with channel indices in case of different sampling frequencies (only for EDF)
%   'chantype'       = string or cell-array with strings, channel types to be read (only for NeuroOmega and BlackRock)
%   'coordsys'       = string, 'head' or 'dewar' (default = 'head')
%   'headerformat'   = name of a MATLAB function that takes the filename as input (default is automatic)
%   'password'       = password structure for encrypted data set (only for mayo_mef30 and mayo_mef21)
%   'readbids'       = string, 'yes', no', or 'ifmakessense', whether to read information from the BIDS sidecar files (default = 'ifmakessense')
%
% This returns a header structure with the following fields
%   hdr.Fs          = sampling frequency
%   hdr.nChans      = number of channels
%   hdr.nSamples    = number of samples per trial
%   hdr.nSamplesPre = number of pre-trigger samples in each trial
%   hdr.nTrials     = number of trials
%   hdr.label       = Nx1 cell-array with the label of each channel
%   hdr.chantype    = Nx1 cell-array with the channel type, see FT_CHANTYPE
%   hdr.chanunit    = Nx1 cell-array with the physical units, see FT_CHANUNIT
%
% For continuously recorded data, this will return nSamplesPre=0 and nTrials=1.
%
% For some data formats that are recorded on animal electrophysiology
% systems (e.g. Neuralynx, Plexon), the following optional fields are
% returned, which allows for relating the timing of spike and LFP data
%   hdr.FirstTimeStamp      number, represented as 32-bit or 64-bit unsigned integer
%   hdr.TimeStampPerSample  number, represented in double precision
%
% Depending on the file format, additional header information can be
% returned in the hdr.orig subfield.
%
% To use an external reading function, you can specify an external function as the
% 'headerformat' option. This function should take the filename as input argument.
% Please check the code of this function for details, and search for BIDS_TSV as
% example.
%
% The following MEG dataformats are supported
%   CTF (*.ds, *.res4, *.meg4)
%   Neuromag/Elekta/Megin (*.fif)
%   BTi/4D (*.m4d, *.pdf, *.xyz)
%   Yokogawa/Ricoh (*.ave, *.con, *.raw)
%   NetMEG (*.nc)
%   ITAB - Chieti (*.mhd)
%   Tristan Babysquid (*.fif)
%   York Instruments (*.meghdf5)
%
% The following EEG dataformats are supported
%   ANT - Advanced Neuro Technology, EEProbe (*.avr, *.eeg, *.cnt)
%   BCI2000 (*.dat)
%   Biosemi (*.bdf)
%   BrainVision (*.eeg, *.seg, *.dat, *.vhdr, *.vmrk)
%   CED - Cambridge Electronic Design (*.smr)
%   EGI - Electrical Geodesics, Inc. (*.egis, *.ave, *.gave, *.ses, *.raw, *.sbin, *.mff)
%   GTec (*.mat, *.hdf5)
%   Generic data formats (*.edf, *.gdf)
%   Megis/BESA (*.avr, *.swf, *.besa)
%   NeuroScan (*.eeg, *.cnt, *.avg)
%   Nexstim (*.nxe)
%   TMSi (*.Poly5)
%   Mega Neurone (directory)
%   Natus/Nicolet/Nervus (.e files)
%   Nihon Kohden (*.m00, *.EEG)
%   Bitalino OpenSignals (*.txt)
%   OpenBCI (*.txt)
%
% The following spike and LFP dataformats are supported
%   Neuralynx (*.ncs, *.nse, *.nts, *.nev, *.nrd, *.dma, *.log)
%   Plextor (*.nex, *.plx, *.ddt)
%   CED - Cambridge Electronic Design (*.smr)
%   MPI - Max Planck Institute (*.dap)
%   Neurosim  (neurosim_spikes, neurosim_signals, neurosim_ds)
%   Windaq (*.wdq)
%   NeuroOmega (*.mat transformed from *.mpx)
%   Neurodata Without Borders: Neurophysiology (*.nwb)
%
% The following NIRS dataformats are supported
%   Artinis - Artinis Medical Systems B.V. (*.oxy3, *.oxy4, *.oxy5, *.oxyproj)
%   BUCN - Birkbeck college, London (*.txt)
%   SNIRF - Society for functional near-infrared spectroscopy (*.snirf)
%
% The following Eyetracker dataformats are supported
%   EyeLink - SR Research (*.asc)
%   SensoMotoric Instruments - (*.txt)
%   Tobii - (*.tsv)
%
% See also FT_READ_DATA, FT_READ_EVENT, FT_WRITE_DATA, FT_WRITE_EVENT,
% FT_CHANTYPE, FT_CHANUNIT

% Copyright (C) 2003-2021 Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

% TODO channel renaming should be made a general option (see bham_bdf)

persistent cacheheader        % for caching the full header
persistent cachechunk         % for caching the res4 chunk when doing realtime analysis on the CTF scanner
persistent db_blob            % for fcdc_mysql

if isempty(db_blob)
  db_blob = false;
end

if iscell(filename)
  % use recursion to read the header from multiple files
  ft_warning('concatenating header from %d files', numel(filename));
  
  hdr = cell(size(filename));
  for i=1:numel(filename)
    hdr{i} = ft_read_header_INTAN(filename{i}, varargin{:});
  end
  
  allhdr = cat(1, hdr{:});
  if numel(unique([allhdr.label]))==sum([allhdr.nChans])
    % each file has different channels, concatenate along the channel dimension
    for i=1:numel(filename)
      assert(isequal(hdr{i}.Fs, hdr{1}.Fs), 'sampling rates are not consistent over files');
      assert(isequal(hdr{i}.nSamples, hdr{1}.nSamples), 'number of samples is not consistent over files');
      assert(isequal(hdr{i}.nTrials, hdr{1}.nTrials), 'number of trials is not consistent over files');
    end
    combined          = hdr{1}; % copy the first header as the general one
    combined.label    = [allhdr.label];
    combined.chanunit = [allhdr.chanunit];
    combined.chantype = [allhdr.chantype];
    combined.nChans   = sum([allhdr.nChans]);
    combined.orig     = hdr;    % store the original header details of all files
  else
    % each file has the same channels, concatenate along the time dimension
    ntrl = nan(size(filename));
    nsmp = nan(size(filename));
    for i=1:numel(filename)
      assert(isequal(hdr{i}.Fs, hdr{1}.Fs), 'sampling rates are not consistent over files');
      assert(isequal(hdr{i}.label, hdr{1}.label), 'channels are not consistent over files');
      ntrl(i) = hdr{i}.nTrials;
      nsmp(i) = hdr{i}.nSamples;
    end
    % the subsequent code concatenates the files over time, i.e. each file has the same channels
    combined      = hdr{1}; % copy the first header as the general one
    combined.orig = hdr;    % store the original header details of all files
    if all(ntrl==1)
      % each file is a continuous recording
      combined.nTrials  = ntrl(1);
      combined.nSamples = sum(nsmp);
    elseif all(nsmp==nsmp(1))
      % each file holds segments of the same length
      combined.nTrials  = sum(ntrl);
      combined.nSamples = nsmp(1);
    else
      ft_error('cannot concatenate files');
    end
  end % concatenate over channels or over time
  % return the header of the concatenated datafiles
  hdr = combined;
  return
end

% get the options
headerformat   = ft_getopt(varargin, 'headerformat');
retry          = ft_getopt(varargin, 'retry', false);     % the default is not to retry reading the header
chanindx       = ft_getopt(varargin, 'chanindx');         % this is used for EDF with different sampling rates
coordsys       = ft_getopt(varargin, 'coordsys', 'head'); % this is used for ctf and neuromag_mne, it can be head or dewar
coilaccuracy   = ft_getopt(varargin, 'coilaccuracy');     % empty, or a number between 0-2
chantype       = ft_getopt(varargin, 'chantype', {});
password       = ft_getopt(varargin, 'password', struct([]));
readbids       = ft_getopt(varargin, 'readbids', 'ifmakessense');

% this should be a cell array
if ~iscell(chantype); chantype = {chantype}; end

% % optionally get the data from the URL and make a temporary local copy
% filename = fetch_url(filename);

if isempty(headerformat)
  % only do the autodetection if the format was not specified
  headerformat = ft_filetype(filename);
end

if iscell(headerformat)
  % this happens for datasets specified as cell-array for concatenation
  headerformat = headerformat{1};
end

if strcmp(headerformat, 'compressed')
  % the file is compressed, unzip on the fly
  inflated     = true;
  filename     = inflate_file(filename);
  headerformat = ft_filetype(filename);
else
  inflated     = false;
end

% for backward compatibility with https://github.com/fieldtrip/fieldtrip/issues/1585
if islogical(readbids)
  % it should be either yes/no/ifmakessense
  if readbids
    readbids = 'yes';
  else
    readbids = 'no';
  end
end

realtime = any(strcmp(headerformat, {'fcdc_buffer', 'ctf_shm', 'fcdc_mysql'}));

% The checkUniqueLabels flag is used for the realtime buffer in case
% it contains fMRI data. It prevents 1000000 voxel names to be checked
% for uniqueness. fMRI users will probably never use channel names
% for anything.

if realtime
  % skip the rest of the initial checks to increase the speed for realtime operation
  
  checkUniqueLabels = false;
  % the cache and fallback option should always be false for realtime processing
  cache    = false;
  fallback = false;
  
else
  % check whether the file or directory exists
  if  ~exist(filename, 'file')
    ft_error('file or directory ''%s'' does not exist', filename);
  end
  
  checkUniqueLabels = true;
  % get the rest of the options, this is skipped for realtime operation
  cache          = ft_getopt(varargin, 'cache');
  fallback       = ft_getopt(varargin, 'fallback');
  checkmaxfilter = ft_getopt(varargin, 'checkmaxfilter', true);
  
  if isempty(cache)
    if any(strcmp(headerformat, {'bci2000_dat', 'eyelink_asc', 'gtec_mat', 'gtec_hdf5', 'mega_neurone', 'nihonkohden_m00', 'smi_txt', 'biosig'}))
      cache = true;
    else
      cache = false;
    end
  end
  
  % ensure that the headerfile and datafile are defined, which are sometimes different than the name of the dataset
%   cfg.filename     = filename;
%   cfg.headerformat = headerformat;
  [filename, headerfile, datafile] = dataset2files(filename, headerformat);
  if ~strcmp(filename, headerfile) && ~ft_filetype(filename, 'ctf_ds') && ~ft_filetype(filename, 'fcdc_buffer_offline') && ~ft_filetype(filename, 'fcdc_matbin')
    filename     = headerfile;                % this function should read the headerfile, not the dataset
    headerformat = ft_filetype(filename);     % update the filetype
  end
end % if skip initial check

% implement the caching in a data-format independent way
if cache && exist(headerfile, 'file') && ~isempty(cacheheader)
  % try to get the header from cache
  details = dir(headerfile);
  if isequal(details, cacheheader.details)
    % the header file has not been updated, fetch it from the cache
    % fprintf('got header from cache\n');
    hdr = rmfield(cacheheader, 'details');
    
    switch ft_filetype(datafile)
      case {'ctf_ds' 'ctf_meg4' 'ctf_old' 'read_ctf_res4'}
        % for realtime analysis end-of-file-chasing the res4 does not correctly
        % estimate the number of samples, so we compute it on the fly
        sz = 0;
        files = dir([filename '/*.*meg4']);
        for j=1:numel(files)
          sz = sz + files(j).bytes;
        end
        hdr.nTrials = floor((sz - 8) / (hdr.nChans*4) / hdr.nSamples);
    end
    
    return
  end % if the details correspond
end % if cache

% the support for head/dewar coordinates is still limited
if strcmp(coordsys, 'dewar') && ~any(strcmp(headerformat, {'fcdc_buffer', 'ctf_ds', 'ctf_meg4', 'ctf_res4', 'neuromag_fif', 'neuromag_mne'}))
  ft_error('dewar coordinates are not supported for %s', headerformat);
end

% deal with data that is organized according to BIDS
if strcmp(readbids, 'yes') || strcmp(readbids, 'ifmakessense')
  [p, f, x] = fileparts(filename);
  % check whether it is a BIDS dataset with json and tsv sidecar files
  % data in a BIDS tsv file (like physio and stim) will be explicitly dealt with in BIDS_TSV
  isbids = startsWith(f, 'sub-') && ~strcmp(x, '.tsv');
  if isbids
    % try to read the metadata from the BIDS sidecar files
    sidecar = bids_sidecar(filename);
    if ~isempty(sidecar)
      data_json = ft_read_json(sidecar);
    end
    sidecar = bids_sidecar(filename, 'channels');
    if ~isempty(sidecar)
      channels_tsv = ft_read_tsv(sidecar);
    end
    sidecar = bids_sidecar(filename, 'electrodes');
    if ~isempty(sidecar)
      electrodes_tsv = ft_read_tsv(sidecar);
    end
    sidecar = bids_sidecar(filename, 'optodes');
    if ~isempty(sidecar)
      optodes_tsv = ft_read_tsv(sidecar);
    end
    sidecar = bids_sidecar(filename, 'coordsystem');
    if ~isempty(sidecar)
      coordsystem_json = ft_read_json(sidecar);
    end
  end
end

% start with an empty header
hdr = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read the data with the low-level reading function
% please maintain this list in alphabetical order
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch headerformat
   
  case 'nwb'
    ft_hastoolbox('MatNWB', 1);	% when I run this locally outside of ft_read_header it does not work for me
    try
      c = load('namespaces/core.mat');
      nwb_version = c.version;
      nwb_fileversion = util.getSchemaVersion(filename);
      if ~strcmp(nwb_version, nwb_fileversion)
        warning(['Installed NWB:N schema version (' nwb_version ') does not match the file''s schema (' nwb_fileversion{1} '). This might result in an error. If so, try to install the matching schema from here: https://github.com/NeurodataWithoutBorders/nwb-schema/releases'])
      end
    catch
      warning('Something might not be alright with your MatNWB path. Will try anyways.')
    end
    tmp = nwbRead(filename); % is lazy, so should not be too costly
    es_key = tmp.searchFor('ElectricalSeries').keys; % find lfp data, which should be an ElectricalSeries object
%     es_key = es_key(~contains(es_key, 'acquisition'));
    if isempty(es_key)
      error('Dataset does not contain an LFP signal (i.e., no object of the class ''ElectricalSeries''.')
    elseif numel(es_key) > 1 % && isempty(additional_user_input) % TODO: Try to sort this out with the user's help
      % Temporary fix: SpikeEventSeries is a daughter of ElectrialSeries but should not be found here (searchFor update on its way)
      % MOD JESUS 15.09.2022 changed 'lfp' for 'lowpass' to fit INTAN naming
      es_key = es_key(contains(es_key,'lowpass','IgnoreCase',true));
    end
    if numel(es_key) > 1 % in case we weren't able to sort out a single
      error('More than one ElectricalSeries present in data. Please specify which signal to use.')
    else
      eseries = io.resolvePath(tmp, es_key{1});
    end
    if isa(eseries.data, 'types.untyped.DataStub')
      hdr.nSamples = eseries.data.dims(2);
    elseif isa(eseries.data, 'types.untyped.DataPipe')
      hdr.nSamples = eseries.data.internal.maxSize(2);
    else
      warning('Cannot determine number of samples in the data.')
      hdr.nSamples = [];
    end
    hdr.Fs          = eseries.starting_time_rate;
    hdr.nSamplesPre = 0; % for now: hardcoded continuous data
    hdr.nTrials     = 1; % for now: hardcoded continuous data
    hdr.label       = {};
    try
        tmp_ch          = io.resolvePath(tmp, eseries.electrodes.table.path).id.data.load; % electrode names
    catch
        tmp_ch          = io.resolvePath(tmp, eseries.electrodes.path).data.load; % MOD by Jesus 08.09.2022
    end
    for iCh=1:numel(tmp_ch) % TODO: does that work if nwb ids are strings?
      if isnumeric(tmp_ch(iCh))
        hdr.label(iCh,1) = {num2str(tmp_ch(iCh))};
      else
        hdr.label(iCh,1) = tmp_ch(iCh);
      end
    end
    hdr.nChans      = numel(hdr.label);
    [hdr.chanunit{1:hdr.nChans,1}] = deal(eseries.data_unit);
    hdr.chanunit    = strrep(hdr.chanunit, 'volt', 'V');
    hdr.chanunit    = strrep(hdr.chanunit, 'micro', 'u');
    % TODO: hdr.FirstTimeStamp
    % TODO: hdr.TimeStampPerSample
    
    % carry over some metadata
    hdr.orig        = [];
    fn = {'general_experimenter', ...
      'general_institution', ...
      'general_keywords', ...
      'general_lab', ...
      'general_notes', ...
      'general_related_publications', ...
      'general_session_id', ...
      'identifier', ...
      'session_description', ...
      'nwb_version', ...
      'help'};
    for iFn = 1:numel(fn)
      if isprop(tmp, fn{iFn}) && ~isempty(tmp.(fn{iFn}))
        hdr.orig.(fn{iFn}) = tmp.(fn{iFn});
      end
    end
    
end % switch headerformat


% Sometimes, the not all labels are correctly filled in by low-level reading functions. See for example bug #1572.
% First, make sure that there are enough (potentially empty) labels:
if numel(hdr.label) < hdr.nChans
  ft_warning('low-level reading function did not supply enough channel labels');
  hdr.label{hdr.nChans} = [];
end
% Now, replace all empty labels with new name:
if any(cellfun(@isempty, hdr.label))
  ft_warning('channel labels should not be empty, creating unique labels');
  hdr.label = fixlabels(hdr.label);
end

if checkUniqueLabels
  if length(hdr.label)~=length(unique(hdr.label))
    % all channels must have unique names
    ft_warning('all channels must have unique labels, creating unique labels');
    megflag = ft_chantype(hdr, 'meg');
    eegflag = ft_chantype(hdr, 'eeg');
    for i=1:hdr.nChans
      sel = find(strcmp(hdr.label{i}, hdr.label));
      if length(sel)>1
        % renaming the first instance is particularly disruptive when the channels are
        % part of standard MEG or EEG channel set, so that should be avoided
        if any(megflag(sel))
          sel = setdiff(sel, sel(find(megflag(sel), 1)));
        elseif any(eegflag(sel))
          sel = setdiff(sel, sel(find(eegflag(sel), 1)));
        else
          sel = sel(2:end);
        end
        for j=1:length(sel)
          % add a number to the original channel name
          hdr.label{sel(j)} = sprintf('%s-%d', hdr.label{sel(j)}, j);
        end
      end
    end
  end
end

% as of November 2011, the header is supposed to include the channel type (see FT_CHANTYPE,
% e.g. meggrad, megref, eeg) and the units of each channel (see FT_CHANUNIT, e.g. uV, fT)

if ~isfield(hdr, 'chantype') && checkUniqueLabels
  % use a helper function which has some built in intelligence
  hdr.chantype = ft_chantype(hdr);
end

if ~isfield(hdr, 'chanunit') && checkUniqueLabels
  % use a helper function which has some built in intelligence
  hdr.chanunit = ft_chanunit(hdr);
end

% ensure that the output grad is according to the latest definition
if isfield(hdr, 'grad')
  hdr.grad = ft_datatype_sens(hdr.grad);
end

% ensure that the output elec is according to the latest definition
if isfield(hdr, 'elec')
  hdr.elec = ft_datatype_sens(hdr.elec);
end

% ensure that the output opto is according to the latest definition
if isfield(hdr, 'opto')
  try
    hdr.opto = ft_datatype_sens(hdr.opto);
  catch
    % the NIRS optode structure is incomplete when reading/converting it from Homer files
    ft_warning('optode structure is not compliant with FT_DATATYPE_SENS');
  end
end

if (strcmp(readbids, 'yes') || strcmp(readbids, 'ifmakessense')) && isbids
  % the BIDS sidecar files overrule the information that is present in the file header itself
  try
    if exist('data_json', 'var')
      hdr.Fs = data_json.SamplingFrequency;
    end
    if exist('channels_tsv', 'var')
      assert(length(channels_tsv.name)  == hdr.nChans, 'number of channels is not consistent with the BIDS channels.tsv');
      assert(length(channels_tsv.type)  == hdr.nChans, 'number of channels is not consistent with the BIDS channels.tsv');
      assert(length(channels_tsv.units) == hdr.nChans, 'number of channels is not consistent with the BIDS channels.tsv');
      hdr.label     = channels_tsv.name;
      hdr.chantype  = channels_tsv.type;
      hdr.chanunit  = channels_tsv.units;
    end
    if exist('electrodes_tsv', 'var')
      hdr.elec         = [];
      hdr.elec.label   = electrodes_tsv.name;
      hdr.elec.elecpos = [electrodes_tsv.x electrodes_tsv.y electrodes_tsv.z];
    end
    if exist('optodes_tsv', 'var')
      hdr.opto         = [];
      hdr.opto.label   = optodes_tsv.name;
      hdr.opto.optopos = [optodes_tsv.x optodes_tsv.y optodes_tsv.z];
    end
    if exist('coordsystem_json', 'var') && ~isempty(coordsystem_json)
      if isfield(hdr, 'grad')
        if strcmp(coordsys, 'dewar')
          % the sensors will be in dewar coordinates, regardless of the coordsystem_json
        else
          % the coordsystem_json overrules the coordinate system 
          hdr.grad.coordsys = coordsystem_json.MEGCoordinateSystem;
        end
        % the units of the grad structure will be correct, they differ between MEG systems, and may depend on coilaccuracy
        % convert them to the units specified in coordsystem_json
        hdr.grad = ft_convert_units(hdr.grad, coordsystem_json.MEGCoordinateUnits);
      end
      if isfield(hdr, 'elec')
        hdr.elec.coordsys = coordsystem_json.EEGCoordinateSystem;
        hdr.elec.unit     = coordsystem_json.EEGCoordinateUnits;
      end
      if isfield(hdr, 'opto')
        hdr.opto.coordsys = coordsystem_json.NIRSCoordinateSystem;
        hdr.opto.unit     = coordsystem_json.NIRSCoordinateUnits;
      end
    end
  catch ME
    if strcmp(readbids, 'yes')
      ft_error(ME.message);
    else
      ft_warning(ME.message);
    end
  end % catch errors
end % if readbids and isbids

% ensure that these are column arrays and that they do not have empty entries
hdr.label = fixlabels(hdr.label);
if isfield(hdr, 'chantype'), hdr.chantype = fixchantype(hdr.chantype); end
if isfield(hdr, 'chanunit'), hdr.chanunit = fixchanunit(hdr.chanunit); end

% ensure that these are double precision and not integers, otherwise
% subsequent computations that depend on these might be messed up
hdr.Fs          = double(hdr.Fs);
hdr.nSamples    = double(hdr.nSamples);
hdr.nSamplesPre = double(hdr.nSamplesPre);
hdr.nTrials     = double(hdr.nTrials);
hdr.nChans      = double(hdr.nChans);

if inflated
  % compressed file has been unzipped on the fly, clean up
  if strcmp(headerformat, 'brainvision_vhdr')
    % don't delete the header file yet, ft_read_data might still need it
    % the files will be cleaned up by ft_read_data
  else
    delete(filename);
  end
end

if cache && exist(headerfile, 'file')
  % put the header in the cache
  cacheheader = hdr;
  % update the header details (including time stamp, size and name)
  cacheheader.details = dir(headerfile);
  % fprintf('added header to cache\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION to determine the file size in bytes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [siz] = filesize(filename)
l = dir(filename);
if l.isdir
  ft_error('"%s" is not a file', filename);
end
siz = l.bytes;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION to determine the file size in bytes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [hdr] = recursive_read_header(filename)
[p, f, x] = fileparts(filename);
ls = dir(filename);
ls = ls(~strcmp({ls.name}, '.'));  % exclude this directory
ls = ls(~strcmp({ls.name}, '..')); % exclude parent directory
for i=1:length(ls)
  % make sure that the directory listing includes the complete path
  ls(i).name = fullfile(filename, ls(i).name);
end
lst = {ls.name};
hdr = cell(size(lst));
sel = zeros(size(lst));
for i=1:length(lst)
  % read the header of each individual file
  try
    thishdr = ft_read_header_INTAN(lst{i});
    if isstruct(thishdr)
      thishdr.filename = lst{i};
    end
  catch
    thishdr = [];
    ft_warning(lasterr);
    fprintf('while reading %s\n\n', lst{i});
  end
  if ~isempty(thishdr)
    hdr{i} = thishdr;
    sel(i) = true;
  else
    sel(i) = false;
  end
end
sel = logical(sel(:));
hdr = hdr(sel);
tmp = {};
for i=1:length(hdr)
  if isstruct(hdr{i})
    tmp = cat(1, tmp, hdr(i));
  elseif iscell(hdr{i})
    tmp = cat(1, tmp, hdr{i}{:});
  end
end
hdr = tmp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION to fill in empty labels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function labels = fixlabels(labels)
if isnumeric(labels)
  % convert the array of numbers into the corresponding strings
  labels = cellfun(@num2str, num2cell(labels), 'UniformOutput', false);
end
for i = find(cellfun(@isempty, {labels{:}}))
  labels{i} = sprintf('%d', i);
end
labels = labels(:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION to fill in empty chantype
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function chantype = fixchantype(chantype)
if isnumeric(chantype)
  % convert the array of numbers into the corresponding strings
  chantype = cellfun(@num2str, num2cell(chantype), 'UniformOutput', false);
end
sel = cellfun(@isempty, chantype) | strcmp(chantype, 'NaN');
chantype(sel) = {'unknown'};
chantype = chantype(:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION to fill in empty chanunit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function chanunit = fixchanunit(chanunit)
if isnumeric(chanunit)
  % convert the array of numbers into the corresponding strings
  chanunit = cellfun(@num2str, num2cell(chanunit), 'UniformOutput', false);
end
sel = cellfun(@isempty, chanunit) | strcmp(chanunit, 'NaN');
chanunit(sel) = {'unknown'};
chanunit = chanunit(:);
