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
if numel(indices) % subject-specific sessions
    aap=aas_setcurrenttask(aap,modulenum,'subject',indices(1));
else
    aap=aas_setcurrenttask(aap,modulenum);
end

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
                    % check whether the input module(s) has/have been skipped
                    skipped = true(1,numel(deps));
                    sourceaap = aap;
                    sourcenumber = inp.sourcenumber;
                    if sourcenumber == -1 % remote source
                        sourceaap=load(inp.aapfilename);
                        sourceaap=sourceaap.aap;
                        % Store these initial settings before any module specific customisation
                        sourceaap.internal.aap_initial=sourceaap;
                        sourceaap.internal.aap_initial.aap_remote.internal.aap_initial=[]; % Prevent recursively expanding storage
                        if ~isfield(sourceaap.acq_details.subjects(1),'subjname'), sourceaap = aa_convert_subjects(sourceaap); end % aa v < 4.5
                        try
                            sourcenumber = aas_getmoduleindexfromtag(sourceaap, inp.sourcestagename);
                        catch
                            aas_log(aap,true,sprintf('Could not find module in remote aap file with stage tag %s',inp.sourcestagename));
                        end
                    end
                    for d = 1:numel(deps)
                        if ~strcmp(deps{d}{1},inp.sourcedomain) ||... % only relevant modules
                                ~exist(aas_doneflag_getpath_bydomain(sourceaap,inp.sourcedomain,deps{d}{2},sourcenumber),'file') % only finished modules (patch for remote only)
                            continue;
                        end
                        fid = fopen(aas_doneflag_getpath_bydomain(sourceaap,inp.sourcedomain,deps{d}{2},sourcenumber),'r');
                        lines = textscan(fid,'%s\n');
                        fclose(fid);
                        skipped(d) = strcmp(lines{1}{1},'skipped');
                    end
                    if all(skipped) && inp.isessential
                        aas_log(aap,false,sprintf('WARNING: No inputs selected for stream %s. --> MODULE %s will be SKIPPED',inp.name, stagename));
                        close_task(aap,tempdirtodelete);
                        aas_writedoneflag(aap,doneflag,'skipped');
                        return
                    end
                    
                    [gotinputs, streamfiles{inpind}]=aas_retrieve_inputs_part1(aap,inp,allinputs,deps);
                    if isempty(setdiff(gotinputs,allinputs)) && inp.isessential % no new inputs found
                        aas_log(aap,true,sprintf('No inputs obtained for stream %s!\n\tModule %s might not have created it.',inp.name,inp.sourcestagename));
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
close_task(aap,tempdirtodelete);
end

function close_task(aap,tempdirtodelete)
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
end