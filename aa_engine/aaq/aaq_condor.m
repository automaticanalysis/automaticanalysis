classdef aaq_condor<aaq
    properties
        maxretries=5;
        
        filestomonitor=[];
        filestomonitor_jobnum=[];
        compiledfile=[];
        condorpath=[];
        retrynum=[];
        jobnotrun=[];
        jobcount=0;
        jobstatus=[];
        
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
            [aap obj.compiledfile]=make_aws_compiled_tool(obj.aap,which('condor_process_jobq'),compdir) 
            
%            % Clearing all existing Condor jobs for this user
%            aas_log(aap,false,'Clearing all existing Condor jobs for this user');
%            cmd='condor_rm -all';
%            aas_shell(cmd);
            
        end
        %% Queue jobs on Condor:
        %  Write small wrapper textfile (condor_q_job)
        %  Queue job
        %  Watch output files
        
        % Run all tasks on the queue, single threaded
        function [obj]=runall(obj,dontcloseexistingworkers,waitforalljobs)
            global aaparallel
            
            % Now run jobs
            njobs=length(obj.jobqueue);
            obj.retrynum(end+1:njobs)=0;
            obj.jobnotrun(end+1:njobs)=true;
            flaggedretry(1:njobs)=false;
            
            errline={};
            
            itnum=0;
            lsevery=20; % To reduce disc load, only run "ls" every 20 loops, to pick up stragglers
            
            while(any(obj.jobnotrun) || not(isempty(obj.filestomonitor)))
                itnum=itnum+1;
                mon_jobs=obj.jobqueue;
                for i=1:njobs
                    mon_jobs(i).jobnotrun=obj.jobnotrun(i);
                    if (not(obj.fatalerrors) && obj.jobnotrun(i))
                        
                        % Find out whether this job is ready to be allocated by
                        % checking dependencies (done_ flags)
                        readytorun=true;
                        gotdoneflag=false(size(obj.jobqueue(i).tobecompletedfirst));
                        for j=1:length(obj.jobqueue(i).tobecompletedfirst)
                            [tbcpth tbcfle tbcext]=fileparts(obj.jobqueue(i).tobecompletedfirst{j});
                            mon_jobs(i).tobecompletedfirst_status{j}='done';
                            % Check for done file
                            if ~exist(obj.jobqueue(i).tobecompletedfirst{j},'file')
                                % I've no idea why, but NFS or Matlab seems to cache the
                                % directory listing sometimes unless you do
                                % this
                                if mod(itnum,lsevery)==0                                    
                                    [s w]=unix(['ls ' tbcpth]);
                                    % Still not there?
                                    if ~exist(obj.jobqueue(i).tobecompletedfirst{j},'file')
                                        readytorun=false;
                                        mon_jobs(i).tobecompletedfirst_status{j}='waiting';
                                    end
                                else
                                    readytorun=false;
                                    mon_jobs(i).tobecompletedfirst_status{j}='waiting';
                                end;
                            else
                                gotdoneflag(j)=true;
                            end
                        end;
                        
                        % Don't look again for flags we've already found
                        obj.jobqueue(i).tobecompletedfirst(gotdoneflag)=[];
                        
                        if (readytorun)
                            obj.jobcount=obj.jobcount+1;
                            job=obj.jobqueue(i);
                            obj.aap.acq_details.root=aas_getstudypath(obj.aap,job.k);
                            obj.condor_q_job(i,job,obj.retrynum(i));
                            obj.jobnotrun(i)=false;
                            flaggedretry(i)=false;
                        end
                        mon_jobs(i).readytorun=readytorun;
                    end
                end
                % Monitor all of the output files
                donemonitoring=false(size(obj.filestomonitor));
                for ftmind=1:length(obj.filestomonitor)
                    state='initialising';
                    
                    % Start with Condor log file for this job
                    jobnum=obj.filestomonitor_jobnum(ftmind);
                    logfid=fopen(obj.filestomonitor(ftmind).log,'r');
                    if (logfid>0)
                        while(not(feof(logfid)))
                            ln=fgetl(logfid);
                            % Skip any codes that have already been
                            % processed
                            for logcode=1:obj.filestomonitor(ftmind).log_lastcodeprocessed
                                while(not(feof(logfid)) && not(strcmp(deblank(ln),'...')))
                                    ln=fgetl(logfid);
                                end
                                ln=fgetl(logfid);
                            end;
                            
                            % Process the next condor status code
                            if length(ln)>3
%                                 aas_log(obj.aap,false,ln,0.5+0.5*obj.aap.gui_controls.colours.running);
                                switch(str2num(ln(1:3)))
                                    case 0
                                        state='submitted';
                                    case 1
                                        state='executing';
                                    case 5
                                        state='terminated';
                                    case 12
                                        state='held';
                                        aas_log(obj.aap,false,sprintf('Job was held, log %s',ln));
                                        pos_open=strfind(ln,'(');
                                        if ~isempty(pos_open)
                                            pos_close=strfind(ln(pos_open:end),'.');
                                            if ~isempty(pos_close)
                                                pause(0.5); % Pause as probably NFS propagation delays causing error
                                                jobnum=ln(pos_open+1:(pos_open+pos_close-2));
                                                cmd=sprintf('condor_release %s',jobnum);
                                                [s w]=aas_shell(cmd);
                                                if s
                                                    obj.fatalerrors=true;
                                                    aas_log(obj.aap,false,sprintf('Failed to release held job, fatal error.\ncommand: %s\nerror:%s',cmd,w));
                                                else
                                                    aas_log(obj.aap,false,sprintf('Released and marked as retry with %s',cmd));
                                                    obj.retrynum(jobnum)=obj.retrynum(jobnum)+1;
                                                    flaggedretry(jobnum)=true;
                                                end;
                                            end;
                                        end;
                                        
                                end
                                % Don't process this code again
                                obj.filestomonitor(ftmind).log_lastcodeprocessed=obj.filestomonitor(ftmind).log_lastcodeprocessed+1;
                                
                            end;
                            
                            % Find the next code
                            while(not(feof(logfid)) && not(strcmp(deblank(ln),'...')))
                                ln=fgetl(logfid);
                            end
                        end
                        fclose(logfid);
                        
                        obj.jobstatus(jobnum).state=obj.filestomonitor(ftmind).state;
                    end
                    if (not(strcmp(state,obj.filestomonitor(ftmind).state)))
                        %
                        if (strcmp(state,'terminated'))
                            % Look at output file
                            fid=fopen(obj.filestomonitor(ftmind).output);
                            while(not(feof(fid)))
                                ln=fgetl(fid);
                                if (ln==-1)
                                    break;
                                end
                                
                                aas_log(obj.aap,false,ln,obj.aap.gui_controls.colours.running);
                            end
                            fclose(fid);
                            
                            
                            
                            % Look at error file
                            fid=fopen(obj.filestomonitor(ftmind).error);
                            while(not(feof(fid)))
                                ln=fgetl(fid);
                                if (ln==-1)
                                    break;
                                end
                                aas_log(obj.aap,false,ln,'Errors');
                                if not(isempty(deblank(ln)))
                                    if strcmp(deblank(ln),'Killed')
                                        aas_log(obj.aap,false,'Job has been killed, but will assume Condor will restart so not marked as fatal error');
                                    else
                                        if ~flaggedretry(jobnum)
                                            obj.retrynum(jobnum)=obj.retrynum(jobnum)+1;
                                            flaggedretry(jobnum)=true;
                                        end;
                                    end;
                                end
                            end
                            fclose(fid);
                            
                            % Report retry or trigger fatal error
                            if flaggedretry(jobnum)
                                if obj.retrynum(jobnum)>obj.maxretries
                                    obj.fatalerrors=true;
                                    aas_log(obj.aap,true,'Fatal error for job');
                                else
                                    aas_log(obj.aap,false,sprintf('Retry %d/%d for job',obj.retrynum(jobnum),obj.maxretries));
                                    obj.jobnotrun(jobnum)=true;
                                end;
                            end;
                            donemonitoring(ftmind)=true;
                            
                            
                        end
                        
                        aas_log(obj.aap,false,sprintf('PARALLEL (condor) %s:  %s',state,obj.filestomonitor(ftmind).name));
                        obj.filestomonitor(ftmind).state=state;
                        
                    end
                    
                    
                    
                end
                
                mon_filestomonitor=obj.filestomonitor;

                if mod(itnum,20)==0
                    % Save status monitor
                    savejson('mon_jobs',mon_jobs,fullfile(obj.aap.internal.aap_initial.acq_details.root,[obj.aap.internal.aap_initial.directory_conventions.analysisid obj.aap.internal.aap_initial.directory_conventions.analysisid_suffix],'aastatus.json'));
                    savejson('mon_filestomonitor',mon_filestomonitor,fullfile(obj.aap.internal.aap_initial.acq_details.root,[obj.aap.internal.aap_initial.directory_conventions.analysisid obj.aap.internal.aap_initial.directory_conventions.analysisid_suffix],'aastatus.json'));
                end;
                
                % Clear out files we've finished monitoring
                obj.filestomonitor(donemonitoring)=[];
                obj.filestomonitor_jobnum(donemonitoring)=[];
                
                if ~waitforalljobs
                    break;
                end;
                % Lets not overload the filesystem
                pause(0.1);
            end
            if waitforalljobs == 1
                obj.emptyqueue;
            end
            
            if (obj.fatalerrors)
                aas_log(obj.aap,false,'PARALLEL (condor): Fatal errors executing jobs. The errors were:','Errors');
                for errind=1:length(errline)
                    aas_log(obj.aap,false,[' ' errline{errind}],'Errors');
                end;
                aas_log(obj.aap,true,'PARALLEL (condor): Terminating because of errors','Errors');
            end
        end
        
        function [obj]=condor_q_job(obj,jobnum,job,retrynum)
            global aaworker
            
            % Write jobs to data directories to preserve them
            [pth nme ext]=fileparts(job.doneflag);
            condorjobpth=fullfile(obj.condorpath,'condorjob');
            obj.aap=aas_makedir(obj.aap,condorjobpth);
            
            aas_log(obj.aap,false,sprintf('Submitting condor job %s',condorjobpth));
            
            % Get root of job name with date time, attempt number and random id
            [pth nme ext]=fileparts(tempname);
            nme=sprintf('attempt%d_%s_%s',1+retrynum,datestr(now,30),nme);
            
            % Save matlab descriptor of job
            jobfn=fullfile(condorjobpth,[nme '_job.mat']);
            save(jobfn,'job');
            
            % Write condor job submission description
            subfn=fullfile(condorjobpth,[nme '_job.sub']);
            fid=fopen(subfn,'w');
            fprintf(fid,'executable=%s\n',obj.compiledfile);
            fprintf(fid,'universe=vanilla\n');
            fprintf(fid,'environment="MCR_CACHE_ROOT=/tmp%s"\n',condorjobpth);
            
            
            % Check on special requirements
            [stagepath stagename]=fileparts(obj.aap.tasklist.main.module(job.k).name);
            try
                specialrequirements=obj.aap.tasksettings.(stagename)(obj.aap.tasklist.main.module(job.k).index).specialrequirements;
                %    specialrequirements={obj.aap.schema.tasksettings.(stagename)(obj.aap.tasklist.main.module(k).index).ATTRIBUTE.specialrequirements};
            catch
                specialrequirements={};
            end
            
            % Set up filenames of logs
            fles=[];
            fles.log=fullfile(condorjobpth,[ nme '_log.txt']);
            fles.output=fullfile(condorjobpth,[nme '_out.txt']);
            fles.error=fullfile(condorjobpth,[nme '_err.txt']);
            fprintf(fid,'Initialdir=%s\n',condorjobpth);
            fprintf(fid,'log=%s\n',[ nme '_log.txt']);
            fprintf(fid,'output=%s\n',[ nme '_out.txt']);
            fprintf(fid,'error=%s\n',[nme '_err.txt']);
            if isfield(specialrequirements,'imagesize')
                fprintf(fid,'ImageSize=%d\n',specialrequirements.imagesize);
            end;
            if isfield(specialrequirements,'jobtype')
                fprintf(fid,'+JobType="%s"\n',specialrequirements.jobtype);
            end;
            
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
            fles.log_lastcodeprocessed=0;
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
