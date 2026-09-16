function [tf, markerFile] = isSyntheticSubject(subjectFolder)
%ISSYNTHETICSUBJECT  True when a subject folder holds synthetic (fixture) data.
%
% PURPOSE:
%   Synthetic data exists to be analysed, which means it looks exactly like
%   real data to every stage of the pipeline - and that is the danger. One
%   fixture subject swept into a group average silently changes a result that
%   nobody will think to question. So the generator leaves a marker file in the
%   subject's folder, and this is what reads it.
%
%   The marker is a file, not a name: a name can be copied, renamed or guessed
%   at, and a rule based on one ("subjects starting with MR") would either miss
%   a renamed fixture or catch a real animal.
%
% USAGE:
%   tf = isSyntheticSubject('D:\TESTSTUDY\data\preprocessing\MRX')
%   [tf, f] = isSyntheticSubject(folder)
%
% INPUT:
%   subjectFolder - path to the subject's folder in any of the data trees.
%
% OUTPUT:
%   tf         - true when the marker is present.
%   markerFile - the marker's path (whether or not it exists).
%
% SEE ALSO:
%   set_default (drops these subjects from 'all' and refuses to mix them),
%   tests/fixture/makeFakeSession (writes the marker).
%
% Last modified 16.09.2026 (Jesus) - new (#F test fixture).

    markerFile = fullfile(subjectFolder, 'SYNTHETIC_DATA.txt');
    tf = isfile(markerFile);
end
