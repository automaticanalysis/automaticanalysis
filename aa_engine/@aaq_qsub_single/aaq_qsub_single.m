classdef aaq_qsub_single<aaq
    properties
        filestomonitor=[];
        jobnotrun = []
    end
    methods
        function [obj]=aaq_qsub_single(aap)
            obj.aap=aap;
        end
        %% Queue jobs on Qsub:
        %  Queue job
        %  Watch output files
        
        % Run all tasks on the queue, single threaded
        function [obj]=runall(obj,dontcloseexistingworkers,waitforalljobs)
            global aaworker
            
            % Check number of jobs & monitored files
            njobs=length(obj.jobqueue);
            
            % We have already submitted some of these jobs
            submittedJobs = 1:length(obj.jobnotrun);
            obj.jobnotrun = true(njobs,1);
            obj.jobnotrun(submittedJobs) = false;
            
            while any(obj.jobnotrun) && ~waitforalljobs
                
                % Lets not overload the filesystem
                pause(10);
                
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
                            obj.aap.acq_details.root=aas_getstudypath(obj.aap,job.k);
                            % Assign an aap to the job!
                            job.aap=obj.aap;
                            % Run the job
                            obj.qsub_q_job(job);
                            obj.jobnotrun(i)=false;
                        end
                    end
                end
                
                % Don't monitor by default
                donemonitoring=true;
                
                % Monitor all of the output files
                donemonitoring=false(size(obj.filestomonitor));
                
                for ftmind=1:length(obj.filestomonitor)
                    
                    % If output exists, check what it is...
                    if exist(obj.filestomonitor(ftmind).name, 'file')
                        JobLog = load(obj.filestomonitor(ftmind).name);
                        
                        % If we don't have any errors...
                        if ~strcmp(JobLog.optout{3}, 'lasterr')
                            % Check the appropriate columns and print what
                            % happened to the job...
                            
                            moduleName = strtok(JobLog.optout{10}, ':');
                            moduleName = moduleName(strfind(moduleName,'aamod_'):(end-8));
                            
                            aas_log(obj.aap,false,...
                                sprintf('Job %s: %s finished', ...
                                obj.filestomonitor(ftmind).name(1:end-11), moduleName), ...
                                obj.aap.gui_controls.colours.running)
                            
                            aas_log(obj.aap,false,...
                                sprintf('Job used %0.4f hours. and %0.9f GB', ...
                                JobLog.optout{2}./(60*60), JobLog.optout{4}./(1024^3)), ...
                                obj.aap.gui_controls.colours.running)
                            
                            if obj.aap.options.qsub.verbose
                                aas_log(obj.aap,false,...
                                    sprintf('%s', ...
                                    JobLog.optout{10}))
                            end
                            
                            % Also save to file with module name attached!                            
                            fid = fopen(fullfile(aaworker.parmpath,'qsub', 'time_estimates.txt'), 'a');
                            fprintf(fid,'%s\n',moduleName);
                            fprintf(fid,'Job used %0.4f hours. and %0.9f GB\n', ...
                                JobLog.optout{2}./(60*60), JobLog.optout{4}./(1024^3));
                            
                            % Job finished, so no need to monitor
                            donemonitoring(ftmind)=true;
                        else
                            % If a job had an error, it is usually fatal...
                            for e = 1:(length(JobLog.optout{4}.stack)-2)
                                cprintf('r', 'Line %d of %s\n', ...
                                    JobLog.optout{4}.stack(e).line, JobLog.optout{4}.stack(e).file )
                            end                       
                            aas_log(obj.aap,true,...
                                sprintf('Job had an error:\n%s\n', ...
                                JobLog.optout{4}.message), ...
                                obj.aap.gui_controls.colours.running)
                        end
                    end
                end
                
                % Clear out files we've finished monitoring
                obj.filestomonitor(donemonitoring)=[];
                
                % Loop if we are still waiting for jobs to finish...
                if waitforalljobs == 1;
                    if isempty(obj.filestomonitor)
                        waitforalljobs = 0;
                    end
                end
            end            
        end
        
        function [obj]=qsub_q_job(obj,job)
            global aaworker
            
            % Let's store all our qsub thingies in one particular directory
            qsubpath=fullfile(aaworker.parmpath,'qsub');
            if (exist(qsubpath,'dir')==0)
                mkdir(qsubpath);
            end
            cd(qsubpath);
            
            % Submit the job using qsubfeval            
            % Check how much memory and time we should assign to the job
            try
                memReq = obj.aap.tasksettings.(job.stagename).qsub.memoryBase * ... % module specific multiplier
                    obj.aap.options.qsub.memoryMult * ... % study specific multiplier
                    (1024^3); % GB
                timReq = obj.aap.tasksettings.(job.stagename).qsub.timeBase * ... % module specific multiplier
                    obj.aap.options.qsub.timeMult * ... % study specific multiplier
                    60*60; % Hours
            catch
                aas_log(obj.aap,false,...
                    sprintf('%s does not contain information about qsub time/memory requirements!', job.stagename), ...
                    [1 0 0])
                memReq = ... % No module specific multiplier
                    obj.aap.options.qsub.memoryMult * ... % study specific multiplier
                    (1024^3); % GB
                timReq = ... % No module specific multiplier
                    obj.aap.options.qsub.timeMult * ... % study specific multiplier
                    60*60; % Hours
            end
            
            % Submit job
            warning off
            jobid = qsubfeval('aa_doprocessing_onetask', obj.aap,job.task,job.k,job.i,job.j, ... % qsubfeval
                'memreq', memReq, ...
                'timreq', timReq, ...
                'diary', 'always');
            warning on
            % State what the assigned number of hours and GB is...
            fprintf('Job %s, assigned %0.4f hours. and %0.9f GB\n\n', ...
                job.stagename, timReq./(60*60), memReq./(1024^3))
            
            % And monitor for files with the job output
            fles.name=[jobid '_output.mat'];
            fles.state='queued';
            if (isempty(obj.filestomonitor))
                obj.filestomonitor=fles;
            else
                obj.filestomonitor(end+1)=fles;
            end
        end
    end
end
