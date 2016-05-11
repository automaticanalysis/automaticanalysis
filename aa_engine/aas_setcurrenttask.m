% This is used to set up particular components of the aap structure that
% change from module to module. It implements the ability to provide
% module-specific parameters (e.g., for branched pipelines) by applying the
% values in extraparameters.aap to the aap structure after setting it up
%
% function [aap]=aas_setcurrenttask(aap,varargin)
%  1. aap=aas_setcurrenttask(aap) - resets aap to initial state as created by user script
%  2. aap=aas_setcurrenttask(aap,k) - sets aap structure (e.g., root path) to correspond to module k
%  3. aap=aas_setcurrenttask(aap,k,'nodefault') - #2 + keeps global SPM defaults
%  4. aap=aas_setcurrenttask(aap,k,'subject',subj) - #2 + sets subject specific selected session

function [aap]=aas_setcurrenttask(aap,varargin)


% Start with the initial pure aap
initinternal=aap.internal;
initaap=aap.internal.aap_initial;
aap=initaap;
aap.internal=initinternal;

if ~isempty(varargin) && isnumeric(varargin{1})
    k = varargin{1};
    % Set SPM defaults appropriately    
    
    if isfield(aap.schema.tasksettings.(aap.tasklist.main.module(k).name)(aap.tasklist.main.module(k).index).ATTRIBUTE,'modality')
        modality = aap.schema.tasksettings.(aap.tasklist.main.module(k).name)(aap.tasklist.main.module(k).index).ATTRIBUTE.modality;
        switch modality
            case 'MRI'
                aap.spm.defaults.modality = 'FMRI';
                sessions = aap.acq_details.sessions;
                % for backwad compatibilty
                if ~isempty(strfind(aap.tasklist.main.module(k).name,'diffusion'))
                    sessions = 'diffusion_sessions';
                end
            case 'DWI'
                aap.spm.defaults.modality = 'FMRI';
                sessions = aap.acq_details.diffusion_sessions;
            case 'MTI'
                aap.spm.defaults.modality = 'FMRI';
                sessions = aap.acq_details.special_sessions;
            case 'MEG'
                aap.spm.defaults.modality = 'EEG';
                sessions = aap.acq_details.meg_sessions;
            otherwise
                aap.spm.defaults.modality = modality;
                sessions = aap.acq_details.sessions;
        end
    else
        aas_log(aap,0,'WARNING:modality is not set; (F)MRI is assumed');
        modality = 'FMRI';
        aap.spm.defaults.modality = modality; % default modality
        sessions = aap.acq_details.sessions;
    end
    
    if nargin<=2 || ~cell_index(varargin,'nodefault')
        global defaults
        defaults=aap.spm.defaults;
    end
    
    stagename=aas_getstagetag(aap,k);
    
    % Stuff the task likes to know
    aap.tasklist.currenttask.epiprefix=aap.tasklist.main.module(k).epiprefix;
    aap.tasklist.currenttask.extraparameters=aap.tasklist.main.module(k).extraparameters;
    aap.tasklist.currenttask.settings=aap.tasksettings.(aap.tasklist.main.module(k).name)(aap.tasklist.main.module(k).index);
    aap.tasklist.currenttask.inputstreams=aap.schema.tasksettings.(aap.tasklist.main.module(k).name)(aap.tasklist.main.module(k).index).inputstreams;
    aap.tasklist.currenttask.outputstreams=aap.schema.tasksettings.(aap.tasklist.main.module(k).name)(aap.tasklist.main.module(k).index).outputstreams;
    aap.tasklist.currenttask.name=stagename;
    aap.tasklist.currenttask.index=aap.tasklist.main.module(k).index;
    aap.tasklist.currenttask.modulenumber=k;
    aap.tasklist.currenttask.domain=aap.schema.tasksettings.(aap.tasklist.main.module(k).name)(aap.tasklist.main.module(k).index).ATTRIBUTE.domain;
    aap.tasklist.currenttask.modality = modality;
    
    % Recursively copy and parameters in extraparameters.aap for this task on
    % top of aap if provided
    if (isfield(aap.tasklist.main.module(k).extraparameters,'aap'))
        aap=aas_copyparameters(aap.tasklist.main.module(k).extraparameters.aap,aap,'aap');
    end
        
    % Check the apparent study root is set appropriately
    aap.acq_details.root=aas_getstudypath(aap,k);
    % ..and for remote filesystem if we're using one
    remotefilesystem=aap.directory_conventions.remotefilesystem;
    if (~strcmp(remotefilesystem,'none'))
        aap.acq_details.(remotefilesystem).root=aas_getstudypath(aap,remotefilesystem,k);
    else
        aap.acq_details.root=aas_getstudypath(aap,k);
    end
    
    % Check subselected sessions
    selected_sessions=aap.acq_details.selected_sessions;
    if ischar(selected_sessions)
        if strcmp(selected_sessions,'*')
            % Wildcard, same as empty
            selected_sessions=1:numel(sessions);
        else
            % Named sessions, parse to get numbers
            sessionnmes = textscan(selected_sessions,'%s'); sessionnmes = sessionnmes{1};
            selected_sessions=[];
            for sessionnme = sessionnmes'
                sessionind = find(strcmp({sessions.name},sessionnme{1}));
                if isempty(sessionind)
                    aas_log(aap,true,sprintf('Unknown session %s specified in selected_sessions field of a branch in the tasklist, sessions were %s',sessionnme{1},sprintf('%s ',sessions.name)));
                end;
                selected_sessions=[selected_sessions sessionind];
            end
        end
        
    end
    
    if cell_index(varargin,'subject')
        subj = varargin(cell_index(varargin,'subject')+1);
        [junk, subjSess] = aas_getN_bydomain(aap,aas_getsesstype(aap),subj{1});
        selected_sessions = intersect(selected_sessions,subjSess);
    end
    
    aap.acq_details.selected_sessions=selected_sessions;
end

end


