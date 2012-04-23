% rc: 15/12/2009
%  substantial changes to add retrieval of data during multi-stage
%  analysis
%  * needs work from Danny to add this functionality to "internal" looping type


function aa_doprocessing_onetask(aap,task,k,i,j)
global aaworker

aaworker.modulestarttime=now;

if aap.options.timelog
    tic
end

try
    aaworker.parmpath;
catch
    aaworker.parmpath=aaworker_getparmpath(aap,0);
end;

% any task specific settings
aap=aas_setcurrenttask(aap,k);

if (~strcmp(aap.directory_conventions.remotefilesystem,'none'))
    [pth nme ext]=fileparts(tempname);
    tempdirtodelete=datestr(now,30);
    [s w]=system('whoami');
    aap.acq_details.root=fullfile('/cn_data',deblank(w),tempdirtodelete,nme,aap.acq_details.root);
else
    tempdirtodelete=[];
end;

studypath=aas_getstudypath(aap);
aap=aas_makedir(aap,studypath);

% allow full path of module to be provided
[stagepath stagename]=fileparts(aap.tasklist.main.module(k).name);
index=aap.tasklist.main.module(k).index;

if (isfield(aap.schema.tasksettings.(stagename)(index).ATTRIBUTE,'mfile_alias'))
    mfile_alias=aap.schema.tasksettings.(stagename)(index).ATTRIBUTE.mfile_alias;
else
    mfile_alias=stagename;
end;

% retrieve description from module
description=aap.schema.tasksettings.(stagename)(index).ATTRIBUTE.desc;

% find out whether this module needs to be executed once per study, subject or session
domain=aap.schema.tasksettings.(stagename)(index).ATTRIBUTE.domain;

%  If multiple repetitions of a module, add 02,03 etc to end of doneflag
doneflagname=aas_doneflag_getname(aap,k);


% if (rand(1)<0.1)
%     fprintf('On the way to pointless trouble\n');
%     aas_log(aap,true,'Pointless trouble');
% end;

% Clear output stream list
aaworker.outputstreams=[];

% Now execute the module, and change the 'done' flags if task='doit'
switch (domain)
    case 'study'
        [doneflag doneflagpath stagetag]=aas_doneflag_getpath(aap,k);
        studypath=aas_getstudypath(aap);
        aas_makedir(aap,studypath);
        save(fullfile(studypath,sprintf('aap_parameters_%s.mat',stagetag)),'aap');
        if (aas_doneflagexists(aap,doneflag))
            if (strcmp(task,'doit'))
                aas_log(aap,0,sprintf('- completed previously: %s',description),aap.gui_controls.colours.completedpreviously);
                % re-write done flag as it clears dependencies
                aas_writedoneflag(aap,doneflag);
            end;
        else
            switch (task)
                case 'checkrequirements'
                    [aap,resp]=aa_feval(mfile_alias,aap,task);
                    if (length(resp)>0)
                        aas_log(aap,0,['\n***WARNING: ' resp]);
                    end;
                case 'doit'
                    tic
                    % before starting current stage, delete done_
                    % flag for stages that are dependencies of it
                    for k0i=1:length(aap.internal.outputstreamdestinations{k}.stream)
                        aas_delete_doneflag(aap,aap.internal.outputstreamdestinations{k}.stream(k0i).destnumber);
                    end;
                    % now run current stage
                    aas_log(aap,0,sprintf('MODULE %s RUNNING: %s',stagename,description),aap.gui_controls.colours.running);
                    % ...fetch inputs
                    if (~isempty(aap.internal.inputstreamsources{k}))
                        
                        allinputs={};
                        for inpind=1:length(aap.internal.inputstreamsources{k}.stream)
                            inp=aap.internal.inputstreamsources{k}.stream(inpind);
                            
                            % check for inputs at study level
                            gotinputs=aas_retrieve_inputs(aap,inp,allinputs);
                            % ... subject level
                            for i=1:length(aap.acq_details.subjects)
                                gotinputs=[gotinputs; aas_retrieve_inputs(aap,inp,allinputs,i)];
                            end;
                            % ... and session level
                            for i=1:length(aap.acq_details.subjects)
                                for j=aap.acq_details.selected_sessions
                                    gotinputs=[gotinputs; aas_retrieve_inputs(aap,inp,allinputs,i,j)];
                                end;
                            end;
                            if (isempty(gotinputs))
                                aas_log(aap,true,sprintf('No inputs obtained for stream %s',inp.name));
                            end;
                            allinputs=[allinputs;gotinputs];
                        end;
                    end;
                    % ...and run
                    aas_log(aap,false,' executing',aap.gui_controls.colours.executing);
                    [aap,resp]=aa_feval(mfile_alias,aap,task);
                    aas_writedoneflag(aap,doneflag);
                    aas_log(aap,0,sprintf('MODULE %s COMPLETED',stagename),aap.gui_controls.colours.completed);
            end;
            
        end;
    case 'subject'
        [doneflag doneflagpath stagetag]=aas_doneflag_getpath(aap,i,k);
        subjpath=aas_getsubjpath(aap,i);
        aas_makedir(aap,subjpath);
        save(fullfile(subjpath,sprintf('aap_parameters_%s.mat',stagetag)),'aap');
        if (aas_doneflagexists(aap,doneflag))
            if (strcmp(task,'doit'))
                aas_log(aap,0,sprintf('- completed previously: %s for %s',description,aas_getsubjname(aap,i)),aap.gui_controls.colours.completedpreviously);
                %                 re-write done flag as it clears dependencies
                aas_writedoneflag(aap,doneflag);
            end;
        else
            switch (task)
                case 'checkrequirements'
                    [aap,resp]=aa_feval(mfile_alias,aap,task,i);
                    if (length(resp)>0)
                        aas_log(aap,0,['\n***WARNING: ' resp]);
                    end;
                case 'doit'
                    tic
                    % before starting current stage, delete done_
                    % flag for stages that are dependencies of it
                    for k0i=1:length(aap.internal.outputstreamdestinations{k}.stream)
                        aas_delete_doneflag(aap,aap.internal.outputstreamdestinations{k}.stream(k0i).destnumber,i);
                    end;
                    % now run current stage
                    aas_log(aap,0,sprintf('MODULE %s RUNNING: %s for %s',stagename,description,aas_getsubjname(aap,i)),aap.gui_controls.colours.running);
                    
                    % ...fetch inputs
                    if (~isempty(aap.internal.inputstreamsources{k}))
                        
                        allinputs={};
                        for inpind=1:length(aap.internal.inputstreamsources{k}.stream)
                            inp=aap.internal.inputstreamsources{k}.stream(inpind);
                            
                            gotinputs={};
                            % ... and at study level
                            gotinputs=[gotinputs;aas_retrieve_inputs(aap,inp,allinputs)];
                            % ... and at subject level
                            gotinputs=[gotinputs;aas_retrieve_inputs(aap,inp,allinputs,i)];
                            % check for inputs from this stream at session level
                            for j=aap.acq_details.selected_sessions
                                gotinputs=[gotinputs; aas_retrieve_inputs(aap,inp,allinputs,i,j)];
                            end;
                            if (isempty(gotinputs))
                                aas_log(aap,true,sprintf('No inputs obtained for stream %s',inp.name));
                            end;
                            allinputs=[allinputs;gotinputs];
                            
                        end;
                    end;
                    
                    
                    % ...and actually run
                    aas_log(aap,false,' executing',aap.gui_controls.colours.executing);
                    [aap,resp]=aa_feval(mfile_alias,aap,task,i);
                    aas_writedoneflag(aap,doneflag);
                    aas_log(aap,0,sprintf('MODULE %s COMPLETED',stagename),aap.gui_controls.colours.completed);
            end;
            
            
        end;
        
    case 'session'
        [doneflag doneflagpath stagetag]=aas_doneflag_getpath(aap,i,j,k);
        sesspath=aas_getsesspath(aap,i,j);
        aas_makedir(aap,sesspath);
        save(fullfile(sesspath,sprintf('aap_parameters_%s.mat',stagetag)),'aap');
        if (aas_doneflagexists(aap,doneflag))
            if (strcmp(task,'doit'))
                aas_log(aap,0,sprintf('- completed previously: %s for %s',description,aas_getsessname(aap,i,j)),aap.gui_controls.colours.completedpreviously);
                % re-write done flag as it clears dependencies
                aas_writedoneflag(aap,doneflag);
            end;
        else
            switch (task)
                case 'checkrequirements'
                    [aap,resp]=aa_feval(mfile_alias,aap,task,i,j);
                    if (length(resp)>0)
                        aas_log(aap,0,['\n***WARNING: ' resp]);
                    end;
                case 'doit'
                    tic
                    % before starting current stage, delete done_
                    % flag for stages that are dependencies of it
                    for k0i=1:length(aap.internal.outputstreamdestinations{k}.stream)
                        aas_delete_doneflag(aap,aap.internal.outputstreamdestinations{k}.stream(k0i).destnumber,i,j);
                    end;
                    % now run current stage
                    aas_log(aap,0,sprintf('MODULE %s RUNNING: %s for %s',stagename,description,aas_getsessname(aap,i,j)),aap.gui_controls.colours.running);
                    
                    % ...fetch inputs
                    if (~isempty(aap.internal.inputstreamsources{k}))
                        allinputs={};
                        for inpind=1:length(aap.internal.inputstreamsources{k}.stream)
                            inp=aap.internal.inputstreamsources{k}.stream(inpind);
                            gotinputs=aas_retrieve_inputs(aap,inp,allinputs,i,j);
                            gotinputs=[gotinputs;aas_retrieve_inputs(aap,inp,allinputs,i)];
                            gotinputs=[gotinputs;aas_retrieve_inputs(aap,inp,allinputs)];
                            if (isempty(gotinputs))
                                aas_log(aap,true,sprintf('No inputs obtained for stream %s',inp.name));
                            end;
                            allinputs=[allinputs;gotinputs];
                        end;
                    end;
                    % ...and actually run
                    aas_log(aap,false,' executing',aap.gui_controls.colours.executing);
                    [aap,resp]=aa_feval(mfile_alias,aap,task,i,j);
                    aas_writedoneflag(aap,doneflag);
                    aas_log(aap,0,sprintf('MODULE %s COMPLETED',stagename),aap.gui_controls.colours.completed);
            end;
        end;
        
    case 'internal'
        % This domain is for parallelising an internal variable; when
        % running in serial mode, need to run 'parallelise' task to do
        % the global bits, then run 'doit' for each index, or without an index to loop
        % over everything;
        switch (task)
            case 'checkrequirements'
                [aap,resp]=aa_feval(mfile_alias,aap,task);
                if (length(resp)>0)
                    aas_log(aap,0,['***WARNING: ' resp]);
                end;
            case 'doit'
                tic
                % this will be the basis of the path
                [doneflag doneflagpath stagetag]=aas_doneflag_getpath(aap,k);
                
                % i describes the parallel part we're going to run
                if ischar(i)
                    doneflag=fullfile(aas_doneflag_getpath(aap,k),[i '.done']); % file in directory
                    desc=i;
                elseif isstruct(i)
                    try c={i.id};
                    catch c=strcat(struct2cell(i),'.'); % Convert structure to cell
                    end
                    doneflag=fullfile(aas_doneflag_getpath(aap,k),[c{:} 'done']); % file in directory
                    desc=[c{:}];
                end
                if (aas_doneflagexists(aap,doneflag))
                    aas_log(aap,0,sprintf('- completed previously: %s for %s',stagename,desc));
                else
                    
                    % now run current stage
                    aas_log(aap,0,sprintf('MODULE %s RUNNING: for %s',stagename,desc));
                    aas_log(aap,false,' executing');
                    [aap,resp]=aa_feval(mfile_alias,aap,task,i);
                    aas_writedoneflag(aap,doneflag);
                    aas_log(aap,0,sprintf('MODULE %s COMPLETED for %s',stagename,desc));
                end;
        end
        
        
        
    otherwise
        aas_log(aap,1,sprintf('Unknown domain %s associated with stage %s',domain,stagename));
end;

% Tidy up by deleting temporary directory created locally
% This could be shifted to a cache manager
if (tempdirtodelete)
    aap.acq_details.root=tempdirtodelete;
    rmdir(aas_getstudypath(aap),'s');
    % If the directory was changed into this path we'll now get a shell
    % error, so cd
    [s w]=aas_shell('pwd');
    if (s)
        cd ~
    end;
end;

if aap.options.timelog
    time_elapsed
end

%%
function checkmem

return




