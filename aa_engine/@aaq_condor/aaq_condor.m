classdef aaq_condor<aaq
    properties
        filestomonitor=[];
        filestomonitor_jobnum=[];
        compiledfile=[];
        condorpath=[];
    end
    methods
        function [obj]=aaq_condor(aap)
            global aaworker
            obj.aap=aap;
            
            obj.condorpath=fullfile(aaworker.parmpath,'condor');
            if (exist(obj.condorpath,'dir')==0)
                mkdir(obj.condorpath);
            end;
            
            % Compile system
            [pth compdir ext]=fileparts(tempname);
            compdir=fullfile(obj.condorpath,compdir);
            mkdir(compdir);
            [aap obj.compiledfile]=make_aws_compiled_tool(obj.aap,'/cn_developer/camneuro/release-beta-0.0/aws/devel/matlab/condor_process_jobq.m',compdir)
            
            % Clearing all existing Condor jobs for this user
            aas_log(aap,false,'Clearing all existing Condor jobs for this user');
            cmd='condor_rm -all';
            aas_shell(cmd);
            
        end
        %% Queue jobs on Condor:
        %  Write small wrapper textfile (condor_q_job)
        %  Queue job
        %  Watch output files
        
        % Run all tasks on the queue, single threaded
        function [obj]=runall(obj,dontcloseexistingworkers,waitforalljobs)
            global aaparallel
            
            
            
            
            % Now run jobs
            obj.filestomonitor=[];
            njobs=length(obj.jobqueue);
            retrynum=zeros(njobs,1);
            maxretries=5;
            fatalerrors=false;
            errline={};
            jobnotrun=true(njobs,1);
            jobcount=0;
            flaggedretry=false(njobs,1);
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
                            jobcount=jobcount+1;
                            job=obj.jobqueue(i);
                            obj.aap.acq_details.root=aas_getstudypath(obj.aap,job.k);
                            [pth execdir ext]=fileparts(tempname);
                            execdir=fullfile(obj.condorpath,['job_' execdir]);
                            job.aap=obj.aap;
                            jobfn=execdir;
                            save(jobfn,'job');
                            obj.condor_q_job(i,jobfn,job);
                            jobnotrun(i)=false;
                            flaggedretry(i)=false;
                        end
                    end
                end
                % Monitor all of the output files
                donemonitoring=false(size(obj.filestomonitor));
                for ftmind=1:length(obj.filestomonitor)
                    jobnum=obj.filestomonitor_jobnum(ftmind);
                    logfid=fopen(obj.filestomonitor(ftmind).log,'r');
                    if (logfid>0)
                        while(not(feof(logfid)))
                            ln=fgetl(logfid);
                            if length(ln)>3
                                switch(str2num(ln(1:3)))
                                    case 0
                                        state='submitted';
                                    case 1
                                        state='executing';
                                    case 5
                                        state='terminated';
                                end
                            end;
                            while(not(feof(logfid)) && not(strcmp(deblank(ln),'...')))
                                ln=fgetl(logfid);
                            end
                        end
                        fclose(logfid);
                    else
                        state='initialising';
                    end
                    if (not(strcmp(state,obj.filestomonitor(ftmind).state)))
                        %
                        if (strcmp(state,'terminated'))
                            fid=fopen(obj.filestomonitor(ftmind).output);
                            while(not(feof(fid)))
                                ln=fgetl(fid);
                                if (ln==-1)
                                    break;
                                end
                                
                                aas_log(obj.aap,false,ln,obj.aap.gui_controls.colours.running);
                            end
                            fid=fopen(obj.filestomonitor(ftmind).error);
                            while(not(feof(fid)))
                                ln=fgetl(fid);
                                if (ln==-1)
                                    break;
                                end
                                aas_log(obj.aap,false,ln,'Errors');
                                errline{end+1}=ln;
                                if not(isempty(deblank(ln)))
                                    if strcmp(deblank(ln),'Killed')
                                        aas_log(obj.aap,false,'Job has been killed, but will assume Condor will restart so not marked as fatal error');
                                    else
                                        if ~flaggedretry(jobnum)
                                            retrynum(jobnum)=retrynum(jobnum)+1;
                                            flaggedretry(jobnum)=true;
                                        end;
                                    end;
                                end
                            end
                            if flaggedretry(jobnum)
                                if retrynum(jobnum)>maxretries
                                    fatalerrors=true;
                                    aas_log(obj.aap,false,'Fatal error for job');
                                else
                                    aas_log(obj.aap,false,sprintf('Retry %d/%d for job',retrynum(jobnum),maxretries));
                                    jobnotrun(jobnum)=true;
                                end;
                            end;
                            donemonitoring(ftmind)=true;
                        end
                        
                        aas_log(obj.aap,false,sprintf('PARALLEL (condor) %s:  %s',state,obj.filestomonitor(ftmind).name));
                        obj.filestomonitor(ftmind).state=state;
                        
                    end
                    
                    
                    
                end
                
                % Clear out files we've finished monitoring
                obj.filestomonitor(donemonitoring)=[];
                
                % Lets not overload the filesystem
                pause(0.5);
            end
            if waitforalljobs == 1
                obj.emptyqueue;
            end
            
            if (fatalerrors)
                aas_log(obj.aap,false,'PARALLEL (condor): Fatal errors executing jobs. The errors were:','Errors');
                for errind=1:length(errline)
                    aas_log(obj.aap,false,[' ' errline{errind}],'Errors');
                end;
                aas_log(obj.aap,true,'PARALLEL (condor): Terminating because of errors','Errors');
            end
        end
        
        function [obj]=condor_q_job(obj,jobnum,jobfn,job)
            global aaworker
            subfn=tempname;
            fid=fopen(subfn,'w');
            fprintf(fid,'executable=%s\n',obj.compiledfile);
            fprintf(fid,'universe=vanilla\n');
            [pth nme ext]=fileparts(subfn);
            fles=[];
            fles.log=fullfile(obj.condorpath,['log_' nme '.txt']);
            fles.output=fullfile(obj.condorpath,['out_' nme '.txt']);
            fles.error=fullfile(obj.condorpath,['err_' nme '.txt']);
            fprintf(fid,'log=%s\n',fles.log);
            fprintf(fid,'output=%s\n',fles.output);
            fprintf(fid,'error=%s\n',fles.error);
            %            fprintf(fid,'ImageSize=2000000\n');
            fprintf(fid,['arguments="%s %s"\n'],getenv('MCRROOT'), jobfn);
            fprintf(fid,'queue\n');
            fclose(fid);
            % Need to get rid of Matlab libraries from the path, or condor
            % flunks with incompatible libraries
            cmd=sprintf('export LD_LIBRARY_PATH= ;condor_submit %s',subfn);
            [s w]=system(cmd);
            if (s)
                fprintf('Error during condor_submit of %s\n',w);
            end
            
            fles.name=job.description;
            fles.state='queued';
            if (isempty(obj.filestomonitor))
                obj.filestomonitor=fles;
            else
                obj.filestomonitor(end+1)=fles;
            end
            obj.filestomonitor_jobnum(length(obj.filestomonitor))=jobnum;
        end
    end
end
