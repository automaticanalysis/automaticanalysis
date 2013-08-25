
function aas_writedoneflag(aap,fn)
global aaworker
global aaparallel

switch(aap.directory_conventions.remotefilesystem)
    case 's3'
        % Find out how many retries and when
        sql=sprintf('select utctime from %s where analysisid=''%s'' and itemName()=''%s''',aaworker.jobtablename,aap.directory_conventions.analysisid,aaworker.jobid);
        [aap jobentry]=sdb_select(aap,sql);
        % Write to the sdb _doneflags table
        attr=[];
        attr.duration=num2str(24*3600*(now-aaworker.modulestarttime));
        try
            utctimelaunches=jobentry.SelectResult.Item.Attribute.Value;
            attr.utctimelaunches=utctimelaunches;
        catch
        end;
        attr.utctime=num2str(utc_time);
        attr.studypath=aas_getstudypath(aap,aap.directory_conventions.remotefilesystem);
        attr.analysisid=aap.directory_conventions.analysisid;
        attr.jobid=aaworker.jobid;
        attr.processkey=aaparallel.processkey;
        attr.bucket=aaworker.bucket;
        if (isfield(aaworker,'outputstreams') && ~isempty(aaworker.outputstreams))
            osd={};
            for sind=1:length(aaworker.outputstreams)
                st=aaworker.outputstreams(sind);
                osd=[osd sprintf('%s|%d|%s|%s|%s',st.name,st.numoutputs,st.desc,st.s3root,st.filename)];
            end;
            attr.outputstreams=osd;
        end;
        cmd='ifconfig  | grep "inet addr:"| grep -v "127.0.0.1" | cut -d: -f2 | awk "{ print $1}"';
        [s w]=aas_shell(cmd);
        if (~s)
            [ipaddress junk]=strtok(w);
        end;
        attr.ipaddress=ipaddress;
        attr
        [attr success]=sdb_put_attributes(aap,aaworker.doneflagtablename,fn,attr);
        if (~success)
            aas_log(aap,true,sprintf('Error writing done flag %s',fn));
        end;
        % Delete entry from the sdb _jobs table
        %  sdb_delete_attributes(aap,aaworker.jobtablename,fn);
        %      >> Now done in aws_process_jobq
        
        % now need to clear this item's attributes from pending jobs with an sdb query and then attribute deletion...
        %         aas_log(aap,false,sprintf('Analysis ID %s',aap.directory_conventions.analysisid));
        %         aas_log(aap,false,sprintf('Analysis ID suffix %s',aap.directory_conventions.analysisid_suffix));
        
        % For each pending job depedendent on this stage...
        sqlish=sprintf('select * from %s where analysisid=''%s'' and tobecompletedfirst=''%s''',aaworker.jobtablename,[aap.directory_conventions.analysisid aap.directory_conventions.analysisid_suffix],fn);
        fprintf('Getting dependencies with %s\n',sqlish);
        [aap depresp]=sdb_select(aap,sqlish);
        njobs=length(depresp.SelectResult.Item);
        fprintf('Now adjusting %d dependencies\n',njobs);
        for depjobind=1:njobs
            depjobname=depresp.SelectResult.Item(depjobind).Name;
            tmpattr=[];
            tmpattr.tobecompletedfirst=fn;
            
            fprintf('Deleting dependencies of %s\n',fn);
            % Make unique version of this job so SQS & SDB definitely match
            itemname=[depjobname '|' sprintf('%08d',round(rand(1)*1e8))];
            % Remove this dependency from this
            sdb_delete_attributes(aap,aaworker.jobtablename,depjobname,tmpattr);
            % Now, if no remaining tbcf, then we're good to go, add to the
            % queue
%             tbcf_ind=strfind({depresp.SelectResult.Item(depjobind).Attribute.Name},'tobecompletedfirst');
%             numtbcf=sum(~cellfun(@isempty,tbcf_ind))-1;
%             depresp_single=[];
%             depresp_single.SelectResult.Item=depresp.SelectResult.Item(depjobind);
            
            % If more than one dependency, then recheck in case another 
            % was deleted at the same time
%            if (numtbcf>0)

            % Retrieve up-to-date copy in case another thread has already
            % queued it while we've been working through our list
            sqlish=sprintf('select * from %s where itemName()=''%s''',aaworker.jobtablename,depjobname)
            [aap depresp_single]=sdb_select(aap,sqlish);
            tbcf_ind=strfind({depresp_single.SelectResult.Item.Attribute.Name},'tobecompletedfirst');
            numtbcf=sum(~cellfun(@isempty,tbcf_ind));
%            end;
            
            % If no dependencies, then queue it
            if (numtbcf==0)
                % Check it hasn't been queued already by another
                % thread
                q_ind=strfind({depresp_single.SelectResult.Item.Attribute.Name},'queued');
                if (all(cellfun(@isempty,q_ind)))
                    %...so add to the queue of jobs that are ready to go
                    tmpattr=[];
                    tmpattr.queued=num2str(now);
                    sdb_put_attributes(aap,aaworker.jobtablename,depjobname,tmpattr);

                    sessiondetails=[aaworker.username '|' getenv('SESSIONID') '|' getenv('SALT1KEY') '|' getenv('SALT3KEY')];
                    encrypted_parts=aaworker.aacc.crypt.tobase64(aaworker.aacc.crypt.encrypt(int8([aap.directory_conventions.analysisid '|' itemname])));
                    jobpath_encrypted=[sessiondetails '|' char(encrypted_parts(:)') ];
%                     jobpath_encrypted=aaworker.aacc.crypt.tobase64(int8([sessiondetails '|' char(encrypted_parts(:)') ]));
%                 sqs_send_message(aap,'cnanon_jobs',jobpath_encrypted);
%    
%                     jobpath_encrypted=aaworker.aap.crypt.b64_encoder.encode([aaworker.username '|' getenv('SESSIONID') '|' getenv('SALT1KEY') '|' getenv('SALT3KEY') '|' aaworker.aacc.crypt.encrypt(int8([aap.directory_conventions.analysisid  '|' itemname])) ]);
                    
                    sqs_send_message(aap,'cnanon_jobs',jobpath_encrypted);
                    aas_log(aap,false,sprintf('Job %s now queued as good-to-go',depjobname));
                else
                    aas_log(aap,false,sprintf('Job %s already queued',depjobname));
                end;
            else
                aas_log(aap,false,sprintf('Job %s is still dependent on other stages',depjobname));
            end;
        end;
    otherwise
        

        [pth nme ext]=fileparts(fn);
        aas_makedir(aap,pth);
        fid=fopen(fn,'w');
        if (~fid)
            aas_log(aap,1,['Error writing done flag ' fn])
        end;
        
        % Execution time
        fprintf(fid,'%f\n',toc);
        
        % And dump IP of this machine to done flag
        cmd='ifconfig  | grep "inet addr:"| grep -v "127.0.0.1" | cut -d: -f2 | awk "{ print $1}"';
        [s w]=aas_shell(cmd);
        if (~s)
            [ipaddress junk]=strtok(w);
            fprintf(fid,'%s\n',ipaddress);
        end;
        fclose(fid);
        
end;
