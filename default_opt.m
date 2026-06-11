function opts = default_opt()
% default_opt  Return the canonical default options struct for the NGL toolbox.
%
% PURPOSE:
%   Schema-driven thin wrapper. The real source of truth is
%   functions/config/optSchema.m — each option is one entry there
%   (name, default, validator, group, help). This file generates the
%   defaults struct on demand and returns it in the same shape every
%   caller has used historically.
%
% USAGE:
%   opts = default_opt();
%   Called internally by set_default. Users should not call this directly;
%   instead, set options in the opt struct inside NGL_SetAndRunMe.m, or
%   consult one option with opt_help('popDyn.pcaConditions').
%
% OUTPUT:
%   opts - struct with every NGL option field pre-filled to its safe
%          default. Field names and shapes are stable for backward
%          compatibility with downstream callers (set_default Section 2
%          validators, NGL02 / NGL04_* scripts, etc.). verifySchemaParity
%          asserts this matches the previous hand-maintained body.
%
% ADDING A NEW OPTION:
%   Append a single optEntry(...) line to functions/config/optSchema.m
%   in the right group. That is the only file you touch when adding an
%   option going forward. The defaults appear here automatically;
%   set_default's per-entry validation (when fully schema-driven) will
%   pick up the validator without further edits.
%
% Last modified 09.06.2026 (Jesus) - schema-driven via optSchema.

    % functions/config is not on the MATLAB path until set_default's
    % Section 7 runs, but default_opt is called from set_default's
    % Section 1. Idempotently add the folder here so the chain works
    % from a cold MATLAB session (i.e. straight after NGL00_Prep).
    here   = fileparts(mfilename('fullpath'));
    cfgDir = fullfile(here, 'functions', 'config');
    if exist(cfgDir, 'dir')
        addpath(cfgDir);
    end

    opts = generateDefaultsFromSchema(optSchema());
end
