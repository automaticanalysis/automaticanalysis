function [aap resp]=aamod_emeg_invertgroup(varargin)
% Do group inversion (seperately for each sensor type if subsequently
% fusing). Parallelised over block type, file stem (epoch type), and
% sensor type.
%
% By applying priors from group level to individual inversions, the
% accuracy of each individal inversion may decline, but the consistency
% across the group should improve. See Litvak & Friston 2008.
%
% Will clear any exisiting inversions. See aamod_emeg_inversion for code to
% check for existing inversions without repeating unnecessarily. A bit
% complicated to tranfer this to conext of group inversions, so safer just
% to start from scratch.
%
% NB: group inversion only makes a difference when multiple priors are used
% e.g. MSP but not IID
%
% Danny Mitchell 02/04/08; 07/11/08
% Adapted from Rik's cbu_meeg_spm5_pipeline.m
%
% To do:
% File collection needs to be generalised to accept arbitrary filters
% Allow inversion of EEG
% Consider inverting events seperately?
% Consider inverting time windows seperately?
% Consider way of dealing with multiple inversion types? multiple forward
% models? and how/whether to overwrite existing inversions
% Allow user choice of hpf/lpf/hanning etc.

%% check task settings, subject, block etc
[aap subblock doit resp settings]=aa_emeg_checktasksettings(mfilename('fullpath'),varargin);
if ~doit; return; end

%% Parallelise %%%%%%%%%%%%%%%% (prepare globals and looping variable here)
gfile=fullfile(aap.acq_details.root,...
    sprintf('%s_%g.parallel.mat',mfilename,(aap.tasklist.currenttask.index)));

if isempty(varargin) || strcmp(varargin{2},'parallelise')

    % define channel types
    switch settings.SensorTypes; 
        case 1, typ={'-mags','-grds'};
        case 2, typ={''};
        case 3, typ={'-mags','-grds',''};
    end
    Ntyp=length(typ);

    % collect file names
    fprintf('\nCollecting files for each subject, block type, epoch length and sensor type:\n')
    for sub=1:length(aap.acq_details.subjects)
        for sess=aap.acq_details.selected_sessions
            for typn=1:Ntyp
                if isempty(regexp(typ{typn},{'-mags','-grds'},'ONCE'))
                    % no 'split' prefix
                    filt=sprintf('^ace.*_BLOCK.*%s\\.mat$',typ{typn});
                else % include split prefix
                    filt=sprintf('^sace.*_BLOCK.*%s\\.mat$',typ{typn});
                end
                files=aas_emeg_findfiles(aap,filt,{sub,sess});
                for fstem=1:length(files)
                    finalnam{sub,sess,fstem,typn}=files{fstem};
                end
            end
        end
        fprintf('.')
    end

    % prepare looping variable
    loopvar={};
    for sess=aap.acq_details.selected_sessions
        for typn=1:Ntyp
            for fstem=1:length(files)
                loopvar{end+1}=struct('subfiles',char(finalnam(:,sess,fstem,typn)), ...
                    'sess',aap.acq_details.sessions(sess).name, ...
                    'fstem',fstem, ...
                    'sensortype',typ{typn}, ...
                    'id',sprintf('%s.%g.%s',aap.acq_details.sessions(sess).name,fstem,typ{typn}));
            end
        end
    end
    save(gfile,'loopvar');
    %aap.tasksettings.(mfilename).PARALLEL=gfile;

    if ~isempty(varargin); return; end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Let's get ready to run
try load(gfile); loopvar;
catch aas_log(aap,1,sprintf('\n***Failed to prepare parallelisation structure!\n%s',resp));
end;

if isempty(subblock); subblock=1:length(loopvar); % all jobs in serial mode
else subblock=subblock{1}; % just run specified job
end

try settings=rmfield(settings,'specialrequirements'); catch; end % other settings will be converted to ID string

for e=subblock
    %% Do group inversion
    %%% for now I'll just allow single inversion type, inverting whole epoch
    %%% and with event types being inverted together
    
    % build ID string for these inversion settings
    temp=settings;
    try temp=rmfield(temp,'specialrequirements'); catch; end
    try temp=rmfield(temp,'Overwrite'); catch; end
    try temp=rmfield(temp,'timeadded'); catch; end
    try temp=rmfield(temp,'ContrastType'); catch; end
    temp=struct2cell(orderfields(temp));
    id='';
    for c=2:length(temp);
        id=[id regexprep(num2str(temp{c}),'\s*',',') '_'];
    end
    id=[regexprep(id(1:end-1),{'/','\'},{'',''})];
                    
    clear DD temp c varargin gfile doit
    fprintf('\nProcessing: %s block, file %g, %s:',loopvar{e}.sess,loopvar{e}.fstem,loopvar{e}.sensortype);
    fprintf('\n - Loading data');
    subfiles=cellstr(loopvar{e}.subfiles);
    for sub=1:length(subfiles)
        % load data
        DD{sub} = spm_eeg_ldata(subfiles{sub});
        DD{sub}.val = 1; % put inversions into first slot of output
        DD{sub}.inv{1}.comment=id;
        fprintf('.')
        
        % try to save some memory! (assumes single letter prefix
        % convention)
        DD{sub}.ica=sprintf('See %s',DD{sub}.fname(2:end));
        DD{sub}.thresholds=sprintf('See %s',DD{sub}.fname(2:end));
        DD{sub}.MaxFilter=sprintf('See %s',DD{sub}.fname(2:end));
        DD{sub}.filter=sprintf('See %s',DD{sub}.fname(2:end));
        DD{sub}.baseline=single(DD{sub}.baseline);
        DD{sub}.scale=single(DD{sub}.scale);
        DD{sub}.events.reject=logical(DD{sub}.events.reject);
        DD{sub}.events.code=int16(DD{sub}.events.code);
        DD{sub}.events.time=int32(DD{sub}.events.time);
        DD{sub}.channels.eeg=uint16(DD{sub}.channels.eeg);
        DD{sub}.channels.order=uint16(DD{sub}.channels.order);
        DD{sub}.channels.Loc=sprintf('See %s',DD{sub}.fname(2:end));
        DD{sub}.channels.Weight=sprintf('See %s',DD{sub}.fname(2:end));
        DD{sub}.channels.scaled=sprintf('See %s',DD{sub}.fname(2:end));
        DD{sub}.channels.Orient=sprintf('See %s',DD{sub}.fname(2:end));
        DD{sub}.channels.posD=sprintf('See %s',DD{sub}.fname(2:end));
        DD{sub}.channels.ort3D=sprintf('See %s',DD{sub}.fname(2:end)); 
        for v=1:length(DD{sub}.inv)
            DD{sub}.inv{v}.datareg=sprintf('See %s',DD{sub}.fname(2:end));
            DD{sub}.inv{v}.forward.OPTIONS=sprintf('See %s',DD{sub}.fname(2:end));
            % this last one is massive!
        end              
    end
    clear sub subfiles
     
    % define inversion parameters
    if settings.RotatingDipoles
        DD{1}.inv{1}.forward.gainmat=DD{sub}.inv{v}.forward.gainxyz;
    end
    DD{1}.inv{1}.inverse = [];
    % DD{1}.inv{1}.inverse.trials= []; % omit this field to invert all event types
    DD{1}.inv{1}.inverse.type  = settings.Method; % e.g. 'GS'
    DD{1}.inv{1}.inverse.smooth= settings.SourceSmoothness;
    DD{1}.inv{1}.inverse.Np    = settings.SparsePriorsPerHemisphere; % See Friston et al submitted 07. Per hemisphere...so total=Np*3? Also depends on mesh resolution?
    DD{1}.inv{1}.inverse.woi   = []; % Invert whole epoch (is there a risk of activity being displaced through time???)
    DD{1}.inv{1}.inverse.lpf   = -Inf;   % Confusing, but hpf refers to cutoff
    DD{1}.inv{1}.inverse.hpf   = 40;
    DD{1}.inv{1}.comment    = id;
    DD{1}.inv{1}.inverse.Nm = [];   % Reset any parameters from last inversion
    DD{1}.inv{1}.inverse.Nr = [];   % Reset any parameters from last inversion
    DD{1}.inv{1}.inverse.Han = 0;
        
    fprintf('\nRunning group inversion:\n');
    DD = spm_eeg_invert(DD); % not yet saved...

    % Save outputs
    fprintf('\nSaving results');
    for sub=1:length(DD)
        [pth nam]=fileparts(DD{sub}.fname);
        DD{sub}.fname=[nam '_grp.mat'];
        outname = fullfile(DD{sub}.path,DD{sub}.fname);
        %grpfinalnam{sub,sess,fstem,typn}=fullfile(D.path, D.fname);
        if ~exist(outname,'file')
            D = DD{sub};
        else
            tries=0;
            while tries<30
                try rehash; load(outname); break
                catch; % have got temporary errors here;
                    % not sure why; might want to wait and try again...
                    pause(2); tries=tries+1;
                end
            end
            % find or append index for this inversion type (will be
            % overwritten)
            status='unknown';
            for v=1:length(D.inv)
                if ~isfield(D.inv{v},'comment') || isempty(D.inv{v}.comment)
                    status='empty'; break; % comment field has not been
                    % created or filled, so overwrite
                end
                if strcmp(D.inv{v}.comment,id);
                    status='found';break
                end
            end
            if strcmp(status,'unknown'); v=v+1; end
            if isempty(v); v=1; end
            D.inv{v}=DD{sub}.inv{1};
            D.inv{v}.comment=id;
        end
        save(outname, 'D');
        fprintf('.')

        %%% could do contrasts here, but see seperate module

    end % end sub

end % next block type; file stem (epoch type); sensor type

return
