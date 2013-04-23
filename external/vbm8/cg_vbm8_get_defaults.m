function varargout = cg_vbm8_get_defaults(defstr, varargin)
% Get/set the defaults values associated with an identifier
% FORMAT defval = cg_vbm8_get_defaults(defstr)
% Return the defaults value associated with identifier "defstr". 
% Currently, this is a '.' subscript reference into the global  
% "defaults" variable defined in spm_defaults.m.
%
% FORMAT cg_vbm8_get_defaults(defstr, defval)
% Sets the vbm8 value associated with identifier "defstr". The new
% vbm8 value applies immediately to:
% * new modules in batch jobs
% * modules in batch jobs that have not been saved yet
% This value will not be saved for future sessions of SPM. To make
% persistent changes, edit cg_vbm8_defaults.m.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% based on Volkmar Glauches version of
% spm_get_defaults
% $Id: cg_vbm8_get_defaults.m 404 2011-04-11 10:03:40Z gaser $

global vbm8;
if isempty(vbm8)
    cg_vbm8_defaults;
end

% construct subscript reference struct from dot delimited tag string
tags = textscan(defstr,'%s', 'delimiter','.');
subs = struct('type','.','subs',tags{1}');

if nargin == 1
    varargout{1} = subsref(vbm8, subs);
else
    vbm8 = subsasgn(vbm8, subs, varargin{1});
end
