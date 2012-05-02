classdef aaq_qsub<aaq
    properties
        filestomonitor=[];
    end
    methods
        function [obj]=aaq_qsub(aap)
            obj.aap=aap;
        end
        %% Queue jobs on Qsub:
        %  Stop previous workers if these exist
        %  Queue job
        %  Watch output files (???)
        
        % Run all tasks on the queue, single threaded
        function [obj]=runall(obj,dontcloseexistingworkers)
            
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
            jobid = qsubfeval('aa_doprocessing_onetask', obj.aap,job.task,job.k,job.i,job.j, ... % qsubfeval
                'memreq', 4*(1024^3), ...
                'timreq', 1*24*60*60, ...
                'diary', 'always');
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
