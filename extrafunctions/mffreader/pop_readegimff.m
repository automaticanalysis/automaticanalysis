% pop_readegimff() - import an EGI Meta File Format (MFF) recording (pop out window if no arguments).
%
% Usage:
%   >> EEG = pop_readegimff;             % a window pops up
%   >> EEG = pop_readegi( filename );
%   >> EEG = pop_readegi( filename, dtype, chanlocfile );
%
% Inputs:
%   filename       - EGI file name
%   dtype          - Type of data to read ('EEG' or 'PIB')
%   chanlocfile        - [string] channel location file name. Default is
%                    'auto' (autodetection)
%
% Outputs:
%   EEG            - EEGLAB data structure
%
% See also:
%   eeglab(), readegimff(), pop_readegi()
%
% Copyright (C) 2011 Srivas Chennu, University of Cambridge, srivas@gmail.com
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

function [EEG, command] = pop_readegimff(filename, varargin)

EEG = [];
command = '';

if nargin < 1
    % ask user
    [filename, filepath] = uigetfile('*','Choose an EGI MFF file -- pop_readegimff()');
    drawnow;
    if filename == 0
        return;
    end;
    filename = [filepath filename];
    
    promptstr    = { 'Data to import (''EEG'' or ''PIB'')' };
    inistr       = { 'EEG' };
    result       = inputdlg2( promptstr, 'Import EGI MFF file -- pop_readegimff()', 1, inistr, 'pop_readegimff');
    if isempty(result); return; end;
    dtype        = result{1};
    if ~strcmp(dtype,'EEG') && ~strcmp(dtype,'PIB')
        error('Type of data to import must be ''EEG'' or ''PIB''');
    end
    
    [head evt data] = readegimff(filename,dtype); % read EGI file header
    
    chanlocfile = '';
    switch head.nchan
        case 129
            chanlocfile = 'GSN-HydroCel-129.sfp';
        case 257
            chanlocfile = 'GSN-HydroCel-257.sfp';
        case 8
            chanlocfile = 'PIB.sfp';
    end;
    promptstr    = { 'Channel location file' };
    inistr       = { chanlocfile };
    result       = inputdlg2( promptstr, 'Import EGI MFF file -- pop_readegimff()', 1, inistr, 'pop_readegimff');
    if isempty(result); return; end;
    chanlocfile      = result{1};
    
else
    param = finputcheck(varargin, {'chanlocfile', 'string', [], 'auto'; ...
        'datatype', 'string', { 'EEG', 'PIB' }, 'EEG'; ...
        'firstepoch', 'integer', [], []; ...
        'lastepoch', 'integer', [], []; ...
        });
    dtype = param.datatype;
    chanlocfile = param.chanlocfile;
    firstepoch = param.firstepoch;
    lastepoch = param.lastepoch;
    
    [head evt data] = readegimff(filename,dtype,firstepoch,lastepoch); % read EGI file header
end

% load data
% ----------
EEG = eeg_emptyset;
EEG.data = data(1:end-1,:,:);

EEG.comments        = [ 'Original file: ' filename ];
EEG.setname 		= 'EGI file';
EEG.nbchan          = head.nchan-1;
EEG.srate           = head.samp_rate;
EEG.trials          = head.segments;
EEG.pnts            = head.segsamps;

EEG = eeg_checkset(EEG);

% importing the events
% --------------------
if ~isempty(evt)
    type ={evt.type}';
    latency = num2cell(cell2mat({evt.sample}')/EEG.srate);
    value = {evt.value}';
    duration = {evt.duration}';
    epoch = num2cell(ceil(cell2mat({evt.sample}')/EEG.pnts));
    
    disp('Importing events.');
    
    if isfield(evt,'codes')
        codes = {evt.codes}';
        EEG = pop_importevent(EEG,'event',cat(2,type,latency,value,duration,epoch,codes),...
            'fields', {'type','latency','value','duration','epoch','codes'},'timeunit',1);
    else
        EEG = pop_importevent(EEG,'event',cat(2,type,latency,value,duration,epoch),...
            'fields', {'type','latency','value','duration','epoch'},'timeunit',1);
    end
    
    EEG = eeg_checkset(EEG, 'makeur');
    EEG = eeg_checkset(EEG, 'eventconsistency');
end

% importing channel locations
% ---------------------------
if ~isempty(chanlocfile)
    disp('Importing channel locations.');
    
    if strcmp(chanlocfile,'auto')
        switch head.nchan
            case 129
                chanlocfile = 'GSN-HydroCel-129.sfp';
            case 257
                chanlocfile = 'GSN-HydroCel-257.sfp';
            case 8
                chanlocfile = 'PIB.sfp';
        end;
    end
    
    chanlocpath = which(chanlocfile);
    if isempty(chanlocpath)
        error('%s not found in path! Make sure the mffreader directory is in your MATLAB path',chanlocfile);
    end
    
    EEG = fixegilocs(EEG,chanlocpath);
end

if nargin < 1
    command = sprintf('EEG = pop_readegi(''%s'', %s);', filename, vararg2str(varargin) );
end;

return;
