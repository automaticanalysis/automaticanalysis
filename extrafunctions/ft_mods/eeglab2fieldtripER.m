function data = eeglab2fieldtripER(EEG,varargin)
% eeglab2fieldtripER() - This function converts epoched EEGLAB datasets to Fieldtrip.
%
% Usage:
%  >> data = eeglab2fieldtripER( EEG, 'key', 'val' ... );
%
% Inputs:
%   EEG         - [struct] EEGLAB structure
%
% Optional inputs (key-value pairs):
%   'reorient'  - detect electrode positions, and reorient e.g. swapped axes (default = 0) (experimental!)
%
% Outputs:
%   data        - FIELDTRIP data structure
%
% Author: Tibor Auer, University of Surrey, June, 2020.

% Copyright (C) 2020 Tibor Auer, University of Surrey, t.auer@surrey.ac.uk
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

cfg = [];
if nargin > 1
    if mod(numel(varargin),2) || ~all(cellfun(@ischar, varargin(1:2:end)))
        ft_error('wrong argumentation');
    end
    cfg = struct(varargin{:});
end
cfg.reorient = ft_getopt(cfg, 'reorient', 0);

EEGLABFILE = fullfile(EEG.filepath,EEG.filename);

% add eeglab ft toolbox if needed
tbpath = '';
lastwarn('')
ft_hastoolbox('eeglab',1);
if ~isempty(strfind(lastwarn,'eeglab'))
    tbpath = strsplit(lastwarn); tbpath = tbpath{2};
end

hdr = ft_read_header(EEGLABFILE);
events = ft_read_event(EEGLABFILE);

if ~isempty(EEG.dipfit)
    data = eeglab2fieldtrip(EEG,'preprocessing', 'dipfit');
else
    data = eeglab2fieldtrip(EEG,'preprocessing', 'none');
end

if ~isfield(data,'trialinfo') % sinlge epoch
    ev = EEG.event([EEG.event.latency] == find(EEG.times==0));
    data.trialinfo = struct2table(keepfields(ev,{'duration','type'}));
    data.trialinfo.type = cellstr(data.trialinfo.type);
    EEG.epoch.event = 1:numel(EEG.event);
    EEG.epoch.eventtype = {EEG.event.type};
    EEG.epoch.eventduration = {EEG.event.duration};
    EEG.epoch.eventlatency = {EEG.event.latency};
    EEG.epoch.eventurevent = {EEG.event.urevent};
end

trlcfg                     = [];
trlcfg.hdr                 = hdr;
trlcfg.trialdef.eventtype  = 'trigger';
trlcfg.trialdef.eventvalue = unique(data.trialinfo.type);
trlcfg.trialdef.prestim    = -data.time{1}(1); % in seconds
trlcfg.trialdef.poststim   = data.time{1}(end)+1/data.fsample; % in seconds
trlcfg.event = events;
trl = round(ft_trialfun_general(trlcfg));
assert(numel(EEG.epoch)==size(data.trialinfo,1),'numel(EEG.epoch) size(data.trialinfo,1) do to match');
if isnumeric(EEG.event(1).type)
    evcompfun = @eq;
elseif ischar(EEG.event(1).type)
    evcompfun = @strcmp;
end
indE = [];
ureventinfo = zeros(0,3);
evs = unique(data.trialinfo.type);
for e = 1:numel(EEG.epoch)
    if numel(EEG.epoch(e).event) == 1
        evmatch = 1;
        if iscell(EEG.epoch(e).eventurevent)
            indUrE = EEG.epoch(e).eventurevent{evmatch};
        else
            indUrE = EEG.epoch(e).eventurevent(evmatch);
        end
        if iscell(EEG.epoch(e).eventurevent)
            currE = EEG.epoch(e).eventtype{evmatch};
        else
            currE = EEG.epoch(e).eventtype;
        end
    else
        evmatch = evcompfun(EEG.epoch(e).eventtype,data.trialinfo.type{e}); evmatch = find(evmatch);
        if sum(evmatch) > 1
            ft_warning('overlap detected in epoch #%d',e)
            evmatch(2:end) = false;
        end
        indUrE = EEG.epoch(e).eventurevent{evmatch};
        currE = EEG.epoch(e).eventtype{evmatch};
    end
    indE = [indE evmatch];
    urevents = EEG.urevent(1:indUrE);
    evIndUnique = sum(strcmp({urevents.type}, currE));
    evIndAll = sum(cellfun(@(x) any(strcmp(evs,x)), {urevents.type}));
    ureventinfo(e,:) = [evIndUnique evIndAll round(urevents(end).latency)/EEG.srate];
end
trl = trl(logical(indE),:);
data = ft_redefinetrial(struct('trl',trl),data);
data.ureventinfo = array2table(ureventinfo,'VariableNames',{'eventnum' 'eventnum_all' 'latency'});

if cfg.reorient
    chind = cellfun(@(x) find([1 strcmp(data.elec.label,x)],1,'last'), {'Fz' 'Cz' 'Oz'})-1;
    if sum(chind~=0) < 2 || ~chind(1)
        ft_warning('Fz and at least one of Cz and Oz must be defined for automatic detection of orientation');
    else
        chind(chind==0) = [];
        EXP_VO = [1 3 2]; % expected order variance -> this should be 1 3 2
        [junk,so] = sort(std(data.elec.elecpos(chind,:))); %
        data.elec.elecpos(:,1:3) = data.elec.elecpos(:,so(EXP_VO));
        data.elec.chanpos(:,1:3) = data.elec.chanpos(:,so(EXP_VO));
    end
end

% clear path
ft_warning('removing %s toolbox from your MATLAB path',tbpath)
rmpath(tbpath)