classdef aaq_qsub<aaq
    properties
        filestomonitor=[];
    end
    methods
        function [obj]=aaq_qsub(aap)
            obj.aap=aap;
        end
        %% Queue jobs on Qsub:
<<<<<<< HEAD
        %  Stop previous workers if these exist
=======
>>>>>>> Improvements to qsub
        %  Queue job
        %  Watch output files (???)
        
        % Run all tasks on the queue, single threaded
        function [obj]=runall(obj,dontcloseexistingworkers)
<<<<<<< HEAD
=======
            global aaworker
>>>>>>> Improvements to qsub
            
            % Check number of jobs & monitored files
            obj.filestomonitor=[];
            njobs=length(obj.jobqueue);
            
            % Check that we have not had fatal errors
            fatalerrors=false;
            
            jobnotrun=true(njobs,1);
            jobcount=0;
            while(any(jobnotrun) || not(isempty(obj.filestomonitor)))
                for i=1:njobs
                    if (not(fatalerrors) && jobnotrun(i))
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
                            jobcount=jobcount+1;
                            job=obj.jobqueue(i);
                            obj.aap.acq_details.root=aas_getstudypath(obj.aap,job.k);
                            % Assign an aap to the job!
                            job.aap=obj.aap;
                            % Run the job
                            obj.qsub_q_job(job);
                            jobnotrun(i)=false;
                        end
                    end
                end
                % Monitor all of the output files
                donemonitoring=false(size(obj.filestomonitor));
                
                for ftmind=1:length(obj.filestomonitor)
                    
                    % If output exists, check what it is...
                    if exist(obj.filestomonitor(ftmind).name, 'file')
                        JobLog = load(obj.filestomonitor(ftmind).name);
<<<<<<< HEAD
                        % Check the appropriate columns and print what
                        % happened to the job...
                        
                        aas_log(obj.aap,false,...
                            sprintf('Job %s finished\n', ...
                            obj.filestomonitor(ftmind).name(1:end-11)), ...
                            obj.aap.gui_controls.colours.running)
                        
                        aas_log(obj.aap,false,...
                            sprintf('Job used %0.0f sec. and %0.9f GB\n', ...
                            JobLog.optout{2}, JobLog.optout{4}./(1024^3)), ...
                            obj.aap.gui_controls.colours.running)
                        
                        % Job finished, so no need to monitor
                        donemonitoring(ftmind)=true;
                        
                        if ~isempty(JobLog.optout{8})
                            % If a job had an error, it is usually fatal...
                            
                            aas_log(obj.aap,false,...
                            sprintf('Job had an error:\n%s\n', ...
                            JobLog.optout{8}), ...
                            obj.aap.gui_controls.colours.running)
=======
                        
                        % If we don't have any errors...
                        if ~strcmp(JobLog.optout{3}, 'lasterr')
                            % Check the appropriate columns and print what
                            % happened to the job...
                            
                            warning off
                            aas_log(obj.aap,false,...
                                sprintf('Job %s finished\n', ...
                                obj.filestomonitor(ftmind).name(1:end-11)), ...
                                obj.aap.gui_controls.colours.running)
                            
                            aas_log(obj.aap,false,...
                                sprintf('Job used %0.4f hours. and %0.9f GB\n', ...
                                JobLog.optout{2}./(60*60), JobLog.optout{4}./(1024^3)), ...
                                obj.aap.gui_controls.colours.running)
                            
                            % Also save to file with module name attached!
                            moduleName = strtok(JobLog.optout{10}, ':');
                            fid = fopen(fullfile(aaworker.parmpath,'qsub', 'time_estimates.txt'), 'a');
                            fprintf(fid,'%s',moduleName(1:end-8));
                            fprintf(fid,'Job used %0.4f hours. and %0.9f GB\n', ...
                                JobLog.optout{2}./(60*60), JobLog.optout{4}./(1024^3));
                            warning on
                            
                            % Job finished, so no need to monitor
                            donemonitoring(ftmind)=true;
                        else
                            % If a job had an error, it is usually fatal...
                            warning off
                            aas_log(obj.aap,false,...
                                sprintf('Job had an error:\n%s\n', ...
                                JobLog.optout{4}.message), ...
                                obj.aap.gui_controls.colours.running)
                            warning on
>>>>>>> Improvements to qsub
                            
                            fatalerrors = true;
                        end
                    end
                end
                
                % Clear out files we've finished monitoring
                obj.filestomonitor(donemonitoring)=[];
                
                % Lets not overload the filesystem
                pause(0.5);
            end
            obj.emptyqueue;
            
            if (fatalerrors)
                aas_log(obj.aap,true,'PARALLEL (qsub): Fatal errors executing jobs');
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
            warning off
<<<<<<< HEAD
            jobid = qsubfeval('aa_doprocessing_onetask', obj.aap,job.task,job.k,job.i,job.j, ... % qsubfeval
                'memreq', 4*(1024^3), ...
                'timreq', 1*24*60*60, ...
                'diary', 'always');
=======
            try
                jobid = qsubfeval('aa_doprocessing_onetask', obj.aap,job.task,job.k,job.i,job.j, ... % qsubfeval
                    'memreq', ...
                    obj.aap.tasksettings.(job.stagename).qsub.memoryBase * ... % module specific multiplier
                    obj.aap.options.qsub.memoryMult * ... % study specific multiplier
                    (1024^3), ... % GB
                    ...
                    'timreq', ...
                    obj.aap.tasksettings.(job.stagename).qsub.timeBase * ... % module specific multiplier
                    obj.aap.options.qsub.timeMult * ... % study specific multiplier
                    60*60, ... % Hours
                    ...
                    'diary', 'always');
            catch
                aas_log(obj.aap,false,...
                    sprintf('%s does not contain information about qsub time/memory requirements!', job.stagename), ...
                    [1 0 0])
                
                jobid = qsubfeval('aa_doprocessing_onetask', obj.aap,job.task,job.k,job.i,job.j, ... % qsubfeval
                    'memreq', ...
                    obj.aap.options.qsub.memoryMult * ... % study specific multiplier
                    (1024^3), ... % GB
                    ...
                    'timreq', ...
                    obj.aap.options.qsub.timeMult * ... % study specific multiplier
                    60*60, ... % Hours
                    ...
                    'diary', 'always');
            end
>>>>>>> Improvements to qsub
            warning on
            
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
