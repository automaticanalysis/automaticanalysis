classdef aaq_qsub<aaq
    properties
        filestomonitor=[];
    end
    methods
        function [obj]=aaq_qsub(aap)
            obj.aap=aap;
        end
        %% Queue jobs on Qsub:
        %  Write small wrapper textfile (qsub_q_job)
        %  Queue job
        %  Watch output files
        
        % Run all tasks on the queue, single threaded
        function [obj]=runall(obj,dontcloseexistingworkers)
            global aaparallel
            
            obj.filestomonitor=[];
            njobs=length(obj.jobqueue);
            
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
                            end;
                        end;
                        
                        if (readytorun)
                            jobcount=jobcount+1;
                            job=obj.jobqueue(i);
                            obj.aap.acq_details.root=aas_getstudypath(obj.aap,job.k);
                            job.aap=obj.aap;
                            jobfn=tempname;
                            save(jobfn,'job');
                            obj.qsub_q_job(jobfn,job);
                            jobnotrun(i)=false;
                        end;
                    end;
                end;
                % Monitor all of the output files
                donemonitoring=false(size(obj.filestomonitor));
                for ftmind=1:length(obj.filestomonitor)
                    
                    logfid=fopen(obj.filestomonitor(ftmind).log,'r');
                    if (logfid>0)
                        while(not(feof(logfid)))
                            ln=fgetl(logfid);
                            switch(str2num(ln(1:3)))
                                case 0
                                    state='submitted';
                                case 1
                                    state='executing';
                                case 5
                                    state='terminated';
                            end;
                            while(not(feof(logfid)) && not(strcmp(deblank(ln),'...')))
                                ln=fgetl(logfid);
                            end;
                        end;
                        fclose(logfid);
                    else
                        state='initialising';
                    end;
                    if (not(strcmp(state,obj.filestomonitor(ftmind).state)))
                        %
                        if (strcmp(state,'terminated'))
                            fid=fopen(obj.filestomonitor(ftmind).output);
                            while(not(feof(fid)))
                                ln=fgetl(fid);
                                if (ln==-1)
                                    break;
                                end;
                                
                                aas_log(obj.aap,false,ln,obj.aap.gui_controls.colours.running);
                            end;
                            fid=fopen(obj.filestomonitor(ftmind).error);
                            while(not(feof(fid)))
                                ln=fgetl(fid);
                                if (ln==-1)
                                    break;
                                end;
                                aas_log(obj.aap,false,ln,'Errors');
                                if (not(isempty(deblank(ln))))
                                    fatalerrors=true;
                                end;
                            end;
                            donemonitoring(ftmind)=true;
                        end;
                        
                        aas_log(obj.aap,false,sprintf('PARALLEL (qsub) %s:  %s',state,obj.filestomonitor(ftmind).name));
                        obj.filestomonitor(ftmind).state=state;
                        
                    end;
                    
                    
                    
                end;
                
                % Clear out files we've finished monitoring
                obj.filestomonitor(donemonitoring)=[];
                
                % Lets not overload the filesystem
                pause(0.5);
            end;
            obj.emptyqueue;
            
            if (fatalerrors)
                aas_log(obj.aap,true,'PARALLEL (qsub): Fatal errors executing jobs');
            end;
        end;
        
        % SUBMIT JOB!!!! @@@@@
        function [obj]=qsub_q_job(obj,jobfn,job)
            global aaworker
            subfn=tempname;
            fid=fopen(subfn,'w');
            [pth nme ext]=fileparts(subfn);
            
            qsubpath=fullfile(aaworker.parmpath,'qsub');
            if (exist(qsubpath,'dir')==0)
                mkdir(qsubpath);
            end;
            fles=[];
            fles.log=fullfile(qsubpath,['log_' nme '.txt']);
            fles.output=fullfile(qsubpath,['out_' nme '.txt']);
            fles.error=fullfile(qsubpath,['err_' nme '.txt']);
            fprintf(fid,'log=%s\n',fles.log);
            fprintf(fid,'output=%s\n',fles.output);
            fprintf(fid,'error=%s\n',fles.error);
            fprintf(fid,'arguments="/usr/local/MATLAB/R2010b %s"\n',jobfn);
            fprintf(fid,'queue\n');
            fclose(fid);
            
            jobid = qsubcellfun(subfn, {1}, ... % qsubfeval
                'memreq', 4*(1024^3), ...
                'timreq', 1*24*60*60, ...
                'diary', 'always');
            
            %{
            if (s)
                fprintf('Error during qsub_submit of %s\n',w);
            end;
            %}
            
            fles.name=job.description;
            fles.state='queued';
            if (isempty(obj.filestomonitor))
                obj.filestomonitor=fles;
            else
                obj.filestomonitor(end+1)=fles;
            end;
        end
    end;
end
