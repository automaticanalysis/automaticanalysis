classdef aaq_qsub<aaq
    properties
        pool = []
        QV = []
    end
    properties (Hidden)
        jobnotrun = []
		taskinqueue = []
        taskstomonitor = []

        % ensure MAXFILTER license
        initialSubmitArguments = '';
    end
    methods
        function [obj]=aaq_qsub(aap)
            global aaworker;
            global aaparallel;
            aaparallel.numberofworkers=1;
            try
                if ~isempty(aap.directory_conventions.poolprofile)
                    profiles = parallel.clusterProfiles;
                    if ~any(strcmp(profiles,aap.directory_conventions.poolprofile))
                        ppfname = which(spm_file(aap.directory_conventions.poolprofile,'ext','.settings'));
                        if isempty(ppfname)
                            aas_log(aap,true,sprintf('ERROR: settings for pool profile %s not found!',aap.directory_conventions.poolprofile));
                        else                            
                            obj.pool=parcluster(parallel.importProfile(ppfname));
                        end
                    else
                        aas_log(aap,false,sprintf('INFO: pool profile %s found',aap.directory_conventions.poolprofile));
                        obj.pool=parcluster(aap.directory_conventions.poolprofile);
                    end
                    switch class(obj.pool)
                        case 'parallel.cluster.Torque'
                            aas_log(aap,false,'INFO: pool Torque is detected');
                            obj.pool.ResourceTemplate = sprintf('-l nodes=^N^,mem=%dGB,walltime=%d:00:00', aaparallel.memory,aaparallel.walltime);
                            if any(strcmp({aap.tasklist.main.module.name},'aamod_meg_maxfilt')) && ... % maxfilt module detected
                                    ~isempty(aap.directory_conventions.neuromagdir) % neuromag specified
                                obj.initialSubmitArguments = ' -W x=\"NODESET:ONEOF:FEATURES:MAXFILTER\"';
                            end
                            obj.pool.SubmitArguments = strcat(obj.pool.SubmitArguments,obj.initialSubmitArguments);
                        case 'parallel.cluster.Generic'
                            aas_log(aap,false,'INFO: Generic engine is detected');
                            obj.pool.IndependentSubmitFcn = obj.SetArg(obj.pool.IndependentSubmitFcn,'walltime',aaparallel.walltime);
                            obj.pool.IndependentSubmitFcn = obj.SetArg(obj.pool.IndependentSubmitFcn,'memory',aaparallel.memory);        
                    end
                else
                    obj.pool = parcluster('local');
                end
                obj.pool.NumWorkers = aaparallel.numberofworkers;
                obj.pool.JobStorageLocation = aaworker.parmpath;
            catch ME
                aas_log(aap,false,'WARNING: Cluster computing is not supported!');
                aas_log(aap,false,sprintf('\tERROR in %s:\n\tline %d: %s',ME.stack(1).file, ME.stack(1).line, ME.message),aap.gui_controls.colours.warning);
                obj.pool=[];
            end
            obj.aap=aap;
        end
        
        function close(obj)
            if ~isempty(obj.pool)
                for j = 1:numel(obj.pool.Jobs)
                    obj.pool.Jobs(j).cancel;
                end
            end
            close@aaq(obj);
        end
        
        %% Queue jobs on Qsub:
        %  Queue job
        %  Watch output files
        
        % Run all tasks on the queue
        function [obj]=runall(obj,dontcloseexistingworkers,waitforalljobs)
            global aaworker
            
            % Check number of jobs & monitored files
            njobs=length(obj.jobqueue);
            
            % We have already submitted some of these jobs
            submittedJobs = 1:length(obj.jobnotrun);
            obj.jobnotrun = true(njobs,1);
            obj.jobnotrun(submittedJobs) = false;
            
            while any(obj.jobnotrun) || waitforalljobs
                
                % Lets not overload the filesystem
                pause(0.1);
                
                for i=1:njobs
                    if (obj.jobnotrun(i))
                        % Find out whether this job is ready to be allocated by
                        % checking dependencies (done_ flags)
                        readytorun=true;
                        for j=1:length(obj.jobqueue(i).tobecompletedfirst)
                            if (~exist(obj.jobqueue(i).tobecompletedfirst{j},'file'))
                                readytorun=false;
                            end
                        end
                        
                        if (readytorun)
                            % Add a job to the queue
                            job=obj.jobqueue(i);
                            % Run the job
                            obj.qsub_q_job(job);
                            obj.jobnotrun(i)=false;
                        end
                    end
                end
                
                taskstarted = [];
                for ftmind=1:numel(obj.taskinqueue)
                    JobID = obj.taskinqueue(ftmind);
                    Jobs = obj.pool.Jobs([obj.pool.Jobs.ID] == JobID);
                    if isempty(Jobs) % cleared by the GUI
                        if obj.QV.isvalid, obj.QV.Hold = false; end
                        obj.fatalerrors = true; % abnormal terminations
                        obj.close;
                        return;
                    end
                    InputArguments = Jobs.Tasks.InputArguments;
                    moduleName = InputArguments{1}.tasklist.main.module(InputArguments{3}).name;
                    state = Jobs.Tasks.State;
                    
                    if ~strcmp(state,'pending')
                        msg = sprintf('MODULE %s RUNNING: Job%d.', moduleName,JobID);
                        aas_log(obj.aap,false,msg);
                        taskstarted(end+1) = ftmind;
                    end
                end
                obj.taskstomonitor = horzcat(obj.taskstomonitor, obj.taskinqueue(taskstarted));
                obj.taskinqueue(taskstarted) = [];
                
                taskreported = [];
                for ftmind=1:numel(obj.taskstomonitor)
                    JobID = obj.taskstomonitor(ftmind);
                    Jobs = obj.pool.Jobs([obj.pool.Jobs.ID] == JobID);
                    if isempty(Jobs) % cleared by the GUI
                        if obj.QV.isvalid, obj.QV.Hold = false; end
                        obj.fatalerrors = true; % abnormal terminations
                        obj.close;
                        return;
                    end
                    Task = Jobs.Tasks;
                    InputArguments = Task.InputArguments;
                    aap = Task.InputArguments{1};
                    moduleName = obj.aap.tasklist.main.module(InputArguments{3}).name;
                    indices = InputArguments{4};
                    datname = ''; datpath = '';
                    if numel(indices) > 0 % subject specified
                        datname = aas_getsubjdesc(aap,indices(1));  
                        datpath = aas_getsubjpath(aap,indices(1)); 
                    end
                    if numel(indices) > 1 % session specified
                        datname = aas_getsessdesc(aap,indices(1),indices(2));  
                        datpath = aas_getsesspath(aap,indices(1),indices(2)); 
                    end
                    
                    state = Task.State;
                    if ~isempty(Task.Error)
                        switch Task.Error.identifier
                            case 'parallel:job:UserCancellation'
                                state = 'cancelled';
                            otherwise
                                state = 'error';
                        end
                    end
                    
                    switch state
                        case 'failed' % failed to launch
                            msg = sprintf('Job%d had failed to launch (Licence?)!\n Check <a href="matlab: open(''%s'')">logfile</a>\n',JobID,...
                                fullfile(obj.pool.JobStorageLocation,Task.Parent.Name,[Task.Name '.log']));
                            % If there is an error, it is fatal...
                            aas_log(obj.aap,true,msg,obj.aap.gui_controls.colours.error)
                            
                            taskreported(end+1) = ftmind;

                        case 'cancelled' % cancelled
                            msg = sprintf('Job%d had been cancelled by user!\n Check <a href="matlab: open(''%s'')">logfile</a>\n',JobID,...
                                fullfile(obj.pool.JobStorageLocation,Task.Parent.Name,[Task.Name '.log']));
                            % If there is an error, it is fatal...
                            aas_log(obj.aap,true,msg,obj.aap.gui_controls.colours.warning)
                            
                            taskreported(end+1) = ftmind;
                        
                        case 'finished' % without error
                            if isempty(Task.FinishTime), continue; end
                            dtvs = dts2dtv(Task.CreateTime);
                            dtvf = dts2dtv(Task.FinishTime);
                            msg = sprintf('JOB %d: \tMODULE %s \tON %s \tSTARTED %s \tFINISHED %s \tUSED %s.',...
                                JobID,moduleName,datname,Task.CreateTime,Task.FinishTime,sec2dts(etime(dtvf,dtvs)));
                            aas_log(obj.aap,false,msg,obj.aap.gui_controls.colours.completed);
                            
                            % Also save to file with module name attached!
                            fid = fopen(fullfile(aaworker.parmpath,'qsub','time_estimates.txt'), 'a');
                            fprintf(fid,'%s\n',msg);
                            fclose(fid);
                            
                            taskreported(end+1) = ftmind;
                            Jobs.delete;
                            
                        case 'error' % running error
                            msg = sprintf('Job%d on <a href="matlab: cd(''%s'')">%s</a> had an error: %s\n',JobID,datpath,datname,Task.ErrorMessage);
                            for e = 1:numel(Task.Error.stack)
                                % Stop tracking to internal
                                if strfind(Task.Error.stack(e).file,'distcomp'), break, end
                                msg = [msg sprintf('<a href="matlab: opentoline(''%s'',%d)">in %s (line %d)</a>\n', ...
                                    Task.Error.stack(e).file, Task.Error.stack(e).line,...
                                    Task.Error.stack(e).file, Task.Error.stack(e).line)];
                            end
                            % If there is an error, it is fatal...
                            obj.fatalerrors = true;
                            aas_log(obj.aap,true,msg,obj.aap.gui_controls.colours.error)
                            
                            taskreported(end+1) = ftmind;
                    end
                end				
                obj.taskstomonitor(taskreported) = [];
                
                % Loop if we are still waiting for jobs to finish...
                if waitforalljobs
                    if isempty(obj.taskinqueue) && isempty(obj.taskstomonitor)
                        waitforalljobs = false;
                    end
                end  
                
                if obj.aap.options.aaworkerGUI
                    % queue viewer
                    % if ~isempty(obj.QV) && ~obj.QV.isvalid % killed
                    %     return
                    % end
                    if ~isempty(obj.pool)
                        if (isempty(obj.QV) || ~obj.QV.OnScreen) % closed
                            obj.QV = aas_qsubViewerClass(obj);
                            obj.QV.Hold = true;
                            obj.QV.setAutoUpdate(false);
                        else
                            obj.QV.UpdateAtRate;
                            if waitforalljobs, obj.QV.Hold = false; end
                        end
                    end
                end
            end
        end
        
        function obj = QVClose(obj)
            if ~isempty(obj.QV) && obj.QV.isvalid
                obj.QV.Close;
                obj.QV.delete;
                obj.QV = [];
            end
        end
        
        function obj = pool_args(obj,varargin)
            global aaparallel;
            memory = aaparallel.memory;
            walltime = aaparallel.walltime;
            
            for iarg = 1:numel(varargin)
                if ~ischar(varargin{iarg}), continue; end
                switch varargin{iarg}
                    case 'mem'
                        if ~isempty(varargin{iarg+1}), memory = varargin{iarg+1}; end
                    case 'walltime'
                        if ~isempty(varargin{iarg+1}), walltime = varargin{iarg+1}; end
                end
            end
            switch class(obj.pool)
                case 'parallel.cluster.Torque'
                    if round(memory) == memory % round
                        memory = sprintf('%dGB',memory);
                    else % non-round --> MB
                        memory = sprintf('%dMB',memory*1000);
                    end
                    obj.pool.SubmitArguments = strcat(sprintf('-q compute -l mem=%s -l walltime=%d',memory,walltime*3600),obj.initialSubmitArguments);
            %                 obj.pool.SubmitArguments = strcat(obj.initialSubmitArguments,...
            %                     sprintf(' -N Mod%02d_',job.k),...
            %                     sprintf('%03d',job.indices));
                case 'parallel.cluster.Generic'
                    aas_log(aap,false,'INFO: Generic engine is detected');
                    obj.pool.IndependentSubmitFcn = obj.SetArg(obj.pool.IndependentSubmitFcn,'walltime',walltime);
                    obj.pool.IndependentSubmitFcn = obj.SetArg(obj.pool.IndependentSubmitFcn,'memory',memory);
            end
        end
        
        function [obj]=qsub_q_job(obj,job)
            global aaworker
            
            % Let's store all our qsub thingies in one particular directory
            qsubpath=fullfile(aaworker.parmpath,'qsub');
            aas_makedir(obj.aap,qsubpath);
            cd(qsubpath);
            
            % Submit the job
            if ~isempty(obj.pool)
                % Check how much memory and time we should assign to the job
                qsubsettings = {'mem',[],'walltime',[]};
%                 if isfield(obj.aap.tasksettings.(job.stagename)(obj.aap.tasklist.main.module(job.k).index),'qsub')
%                     qsub = obj.aap.tasksettings.(job.stagename)(obj.aap.tasklist.main.module(job.k).index).qsub;
%                     for f = fieldnames(qsub)'
%                         switch f{1}
%                             case 'memoryBase'
%                                 qsubsettings{2} = qsub.memoryBase;
%                             case 'timeBase'
%                                 qsubsettings{4} = qsub.timeBase;
%                         end
%                     end
%                 end

                if isa(obj.pool,'parallel.cluster.Torque'), obj = obj.pool_args(qsubsettings{:}); end
                
                J = createJob(obj.pool);
                cj = @aa_doprocessing_onetask;
                nrtn = 0;
                inparg = {obj.aap,job.task,job.k,job.indices, aaworker};
                
                if isprop(J,'AutoAttachFiles'), J.AutoAttachFiles = false; end

                % [RT 2013-09-04 and 2013-11-11; TA 2013-11-14 and 2014-12-12] Make workers self-sufficient by passing
                % them the aa paths. Users don't need to remember to update
                % their own default paths (e.g. for a new aa version)
                global aacache;
                if isprop(J,'AdditionalPaths')
                    J.AdditionalPaths = aacache.path.reqpath;
                elseif isprop(J,'PathDependencies')
                    J.PathDependencies = aacache.path.reqpath;
                end
                
                createTask(J,cj,nrtn,inparg,'CaptureDiary',true);
                J.submit;
                %                 % State what the assigned number of hours and GB is...
                % Naas_movParsot in use [TA]
                %                 fprintf('Job %s, assigned %0.4f hours. and %0.9f GB\n\n', ...
                %                     job.stagename, timReq./(60*60), memReq./(1024^3))

                % And monitor for files with the job output
				obj.taskinqueue(end+1) = J.ID;
            else
                aa_doprocessing_onetask(obj.aap,job.task,job.k,job.indices);
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% UTILS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dtv = dts2dtv(dts)
s = textscan(dts,'%s'); s = s{1}; s(5) = [];
s = strcat(s,{' ',' ',' ',' ',' '}'); s = [s{:}]; s = s(1:end-1);
dtformat = 'ddd mmm dd HH:MM:SS yyyy';
dtv = datevec(s,dtformat);
end

function dts = sec2dts(dt)
dt_str = {'s','m','h'};
dt_div = [60 60 24]; 

dts = '';
for i = 1:numel(dt_str)
    dts = [' ' num2str(mod(dt,dt_div(i))) dt_str{i} dts]; dt = floor(dt/dt_div(i));
    if ~dt, break, end
end
if dt
    dts = [num2str(dt) 'd' dts]; 
else
    dts = dts(2:end);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DLG %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
