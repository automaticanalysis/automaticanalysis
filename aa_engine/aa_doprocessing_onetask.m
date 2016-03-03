% rc: 15/12/2009
%  substantial changes to add retrieval of data during multi-stage
%  analysis
% rc: 3 June 2012
%  changed dependency calculation code to allow for other kinds of parallel
% ta: 13 Nov 2015
%  accept aaworker (qsub)

function aap=aa_doprocessing_onetask(aap,task,modulenum,indices,gaaworker)

global aaworker
if nargin == 5 % aaworker passed
    aaworker = gaaworker;
end

aaworker.modulestarttime=now;
    
if aap.options.timelog
    tic
end

if ~isfield(aaworker,'parmpath')
    aas_log(aap,false,'INFO: No engine detected!\nINFO: aaworker with ID=0 will be created');
    aaworker.parmpath=aaworker_getparmpath(aap,0);
end

% any task specific settings
aap=aas_setcurrenttask(aap,modulenum);

if (~strcmp(aap.directory_conventions.remotefilesystem,'none'))
    [pth nme ext]=fileparts(tempname);
    tempdirtodelete=[datestr(now,30) '_' nme];
    [s w]=system('whoami');
    tempdirtodelete=fullfile('/cn_data',deblank(w),tempdirtodelete);
    aap.acq_details.root=fullfile(tempdirtodelete,aap.acq_details.root);
else
    tempdirtodelete=[];
end;

studypath=aas_getstudypath(aap);
aap=aas_makedir(aap,studypath);

% allow full path of module to be provided
[stagepath stagename]=fileparts(aap.tasklist.main.module(modulenum).name);
index=aap.tasklist.main.module(modulenum).index;

if (isfield(aap.schema.tasksettings.(stagename)(index).ATTRIBUTE,'mfile_alias'))
    mfile_alias=aap.schema.tasksettings.(stagename)(index).ATTRIBUTE.mfile_alias;
else
    mfile_alias=stagename;
end;

taskSchema = aap.schema.tasksettings.(stagename)(index);

% retrieve description from module
description=taskSchema.ATTRIBUTE.desc;

% find out whether this module needs to be executed once per study, subject or session
domain=taskSchema.ATTRIBUTE.domain;

%  If multiple repetitions of a module, add 02,03 etc to end of doneflag
doneflagname=aas_doneflag_getname(aap,modulenum);


% if (rand(1)<0.1)
%     fprintf('On the way to pointless trouble\n');
%     aas_log(aap,true,'Pointless trouble');
% end;

% Clear output stream list
aaworker.outputstreams=[];

% Now execute the module, and change the 'done' flags if task='doit'

[doneflag doneflagpath stagetag]=aas_doneflag_getpath_bydomain(aap,domain,indices,modulenum);
outputpath=aas_getpath_bydomain(aap,domain,indices);
aas_makedir(aap,outputpath);
aap_tmp=aap;
aap.internal.streamcache=[];
save(fullfile(outputpath,sprintf('aap_parameters_%s.mat',stagetag)),'aap');
aap=aap_tmp;
if (aas_doneflagexists(aap,doneflag))
    if (strcmp(task,'doit'))
        aas_log(aap,0,sprintf('- completed previously: %s for %s',description,doneflagpath),aap.gui_controls.colours.completedpreviously);
        %                 re-write done flag as it clears dependencies
        aas_writedoneflag(aap,doneflag);
    end;
else
    switch (task)
        case 'checkrequirements'
            [aap,resp]=aa_feval_withindices(mfile_alias,aap,task,indices);
            
            if (length(resp)>0)
                aas_log(aap,0,['\n***WARNING: ' resp]);
            end;
        case 'doit'
            tic
            % before starting current stage, delete done_
            % flag for stages that are dependencies of it
            for k0i=1:length(aap.internal.outputstreamdestinations{modulenum}.stream)
                aas_delete_doneflag_bydomain(aap,aap.internal.outputstreamdestinations{modulenum}.stream(k0i).destnumber,domain,indices);
            end;
            % now run current stage
            aas_log(aap,0,sprintf('MODULE %s RUNNING: %s for %s',stagename,description,doneflagpath),aap.gui_controls.colours.running);
            
            % ...fetch inputs
            if ~isempty(aap.internal.inputstreamsources{modulenum})
                
                allinputs={};
                streamfiles=cell(length(aap.internal.inputstreamsources{modulenum}.stream),1);
                for inpind=1:length(aap.internal.inputstreamsources{modulenum}.stream)
                    
                    inp=aap.internal.inputstreamsources{modulenum}.stream(inpind);
                    
                    % There might be additional settings for this input 
                    % Added by CW to allow domain override
                    if iscell(taskSchema.inputstreams.stream)  && ...
                            ~(numel(taskSchema.inputstreams.stream) == 1 && isstruct(taskSchema.inputstreams.stream{1}))
                        inputSchema = taskSchema.inputstreams.stream{inpind};
                    else
                        inputSchema = taskSchema.inputstreams.stream{1}(inpind);
                    end

                    % Let's add the ability to force the search domain to
                    % change for a specific input (e.g., for x-val purposes)
                    searchDomain = domain;      % Default to domain of this module
                    searchIndices = indices;    % Default to domain indices
                    if isstruct(inputSchema) && isfield(inputSchema.ATTRIBUTE,'forcedomain')
                        searchDomain = inputSchema.ATTRIBUTE.forcedomain;
                        
                        % The current 'indices' specify the current module,
                        % we have to update them to reflect the new
                        % (forced) search domain.
                        searchDomainTree = aas_dependencytree_finddomain(searchDomain,aap.directory_conventions.parallel_dependencies,{});
                        moduleDomainTree = aas_dependencytree_finddomain(domain,aap.directory_conventions.parallel_dependencies,{});
                        
                        if length(searchDomainTree) < length(moduleDomainTree)
                            searchIndices = searchIndices(1:length(searchDomainTree)-1);
                        else
                            aas_log(aap, 1, 'NYI: forcing domain to be more specific, we have to have a way to specify the new indices.');
                        end
   
                    end
                    
                    deps=aas_getdependencies_bydomain(aap,inp.sourcedomain,searchDomain,searchIndices);
                    % subselect dependency
                    toSelect = logical([]);
                    for d = 1:numel(deps)
                        toSelect(d) = isempty(strfind(deps{d}{1},'session')) ||...      % not a session
                            (numel(deps{d}{2}) < 2) ||...                               % not applicable (e.g. isc_session)
                            any(deps{d}{2}(2) == aap.acq_details.selected_sessions);    % selected session
                    end
                    deps = deps(toSelect);
                    
                    [gotinputs, streamfiles{inpind}]=aas_retrieve_inputs_part1(aap,inp,allinputs,deps);
                    if isempty(gotinputs)
                        aas_log(aap,true,sprintf('No inputs obtained for stream %s',inp.name));
                    end;
                    allinputs=[allinputs;gotinputs];
                end;
                
                % Actually copy the data files
                % This is now done separately to allow for asynchronous retrieval of all of the streams together, where they've  
                % gone of too S3 or Glacier archive
                for inpind=1:length(aap.internal.inputstreamsources{modulenum}.stream)
                    aap=aas_retrieve_inputs_part2(aap,streamfiles{inpind});
                end;
            end;
            
            
            % ...and actually run
            aas_log(aap,false,' executing',aap.gui_controls.colours.executing);
            [aap,resp]=aa_feval_withindices(mfile_alias,aap,task,indices);
            aas_writedoneflag(aap,doneflag);
            aas_log(aap,0,sprintf('MODULE %s COMPLETED',stagename),aap.gui_controls.colours.completed);            
    end;
end;
% Tidy up by deleting temporary directory created locally
% This could be shifted to a cache manager
if (tempdirtodelete)
    rmdir(tempdirtodelete,'s');
    % If the directory was changed into this path we'll now get a shell
    % error, so cd
    [s w]=aas_shell('pwd');
    if (s)
        cd ~
    end;
end;

if aap.options.timelog
    aas_time_elapsed
end

