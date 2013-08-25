% Retrieves inputs from a particular stream to a particular destination
%

function [gotinputs]=aas_retrieve_inputs(aap,inputstream,gotinputs,deps)
global aaworker

for depind=1:length(deps)
    % I'm doing this for you, dear reader. More readable names!
    domain=deps{depind}{1};
    indices=deps{depind}{2};
    
    streamname=inputstream.name;
    sourcestagenumber=inputstream.sourcenumber;
    ismodified=inputstream.ismodified;
    
    % Only take part of stream after last period
    pos=find(streamname=='.');
    if (~isempty(pos))
        fromstreamname=streamname(pos(end)+1:end);
    else
        fromstreamname=streamname;
    end;
    
    %% REMOTE STREAM - this is stored on another machine and retrieved with rsync+ssh
    % Is it a remote stream?
    if (sourcestagenumber==-1)
        % Input data are on a remote machine
        % Get the destination path
        dest=aas_getpath_bydomain(aap,deps{depind}{1},deps{depind}{2});
        
        % Make sure it exists
        [s m mid]=mkdir(dest);
        
        % Don't do this any more - let rsync handle the redundant copies one
        % layer down
        % See if remote files have already been retrieved?
        %remotefetchfn=fullfile(dest,sprintf('done_remote_fetch_of_%s',streamname));
        doneremotefetch=false; %exist(remotefetchfn,'file');
        
        % First we need to get the previous aap file
        remote_aap_fn=fullfile(dest,'remote_aap_parameters.mat');
        if (~doneremotefetch)
            aap=aas_copyfromremote(aap,inputstream.host,inputstream.aapfilename,remote_aap_fn);
        else
            aas_log(aap,false,sprintf('Not repeating earlier fetch of remote data to %s',dest));
        end;
        aap_remote=load(remote_aap_fn);
        aap_remote=aap_remote.aap;
        % Store these initial settings before any module specific customisation
        aap_remote.internal.aap_initial=aap_remote;
        aap_remote.internal.aap_initial.aap_remote.internal.aap_initial=[]; % Prevent recursively expanding storage
        
        % Which module?
        gotit=false;
        for modind=1:length(aap_remote.tasklist.main.module)
            if strcmp(aas_getstagetag(aap_remote,modind),inputstream.sourcestagename)
                gotit=true;
                break
            end;
        end;
        if (~gotit)
            aas_log(aap,true,sprintf('Could not find module in remote aap file with stage tag %s',inputstream.sourcestagename));
        end;
        
        % Set up remote aap for this module
        aap_remote=aas_setcurrenttask(aap_remote,modind);
        
        % Check for recent additions to aap in case we're using an older
        % version on the remote machine
        if ~isfield(aap_remote.directory_conventions,'parallel_dependencies')
            aap_remote.directory_conventions.parallel_dependencies=aap.directory_conventions.parallel_dependencies;
            aas_log(aap,false,'Remote aap does not contain parallel_dependencies field - probably from older aa');
        end;
        % Get remote directory
        remoteindices=indices;
        if length(indices)>2
            findsession=find(strcmp(aap.acq_details.sessions(indices(2)).name,{aap_remote.acq_details.sessions.name}));
            if isempty(findsession)
                aas_log(aap,true,sprintf('Error loading remote file with stage tag %s as could not find session called %s',inputstream.sourcestagename,aap.acq_details.sessions(indices(2)).name));
            end;
            remoteindices(2)=findsession;
        end;
        src=aas_getpath_bydomain(aap_remote,domain,remoteindices,modind);
        
        remoteoutputstreamdesc=fullfile(src,sprintf('stream_%s_outputfrom_%s.txt',fromstreamname,inputstream.sourcestagename));
        
        % This is written by this routine - the input to the stage we're
        % working on
        
        % Delete non-qualified stream name, if this exists, as this will
        % override a qualified filename, which is dangerous
        non_qualified_fn=fullfile(dest,sprintf('stream_%s_inputto_%s.txt',fromstreamname,aap.tasklist.currenttask.name));
        if (exist(non_qualified_fn,'file'))
            delete(non_qualified_fn);
        end;
        
        % Now produce either qualified or non-qualified (depending on
        % streamname)
        inputstreamdesc=fullfile(dest,sprintf('stream_%s_inputto_%s.txt',streamname,aap.tasklist.currenttask.name));
        
        % Copy from one drive to another on local filesystem
        %  Note the src corresponds to the outputstream of the previous
        %  stage, and the dest corresponds to the inputstream of the
        %  next stage. A touch confusing!
        outputstreamdesc=fullfile(dest,sprintf('stream_%s_remoteoutputfrom_%s_%s.txt',fromstreamname,inputstream.host,inputstream.sourcestagename));
        
        if ~doneremotefetch
            aap=aas_copyfromremote(aap,inputstream.host,remoteoutputstreamdesc,outputstreamdesc,true);
        end;
        
        if (exist(outputstreamdesc,'file'))
            
            reloadfiles=true;
            fid=fopen(outputstreamdesc,'r');
            
            % Load checksum with error checking
            [aap md5]=loadmd5(aap,fid,streamname);
            
            % Get filenames
            fns=textscan(fid,'%s');
            fns=fns{1};
            fns_dest=cell(length(fns),1);
            fns_dest_full=cell(length(fns),1);
            wasnamechange=false;
            for ind=1:length(fns)
                % Check to see whether a filename with this name has
                % already been loaded. If so, add unique suffix
                fns_dest{ind}=fns{ind};
                fns_dest_full{ind}=fullfile(dest,fns_dest{ind});
                suffix=1;
                while (1)
                    pos=[strcmp(fns_dest_full{ind},gotinputs)];
                    if (~any(pos))
                        break;
                    end;
                    wasnamechange=true;
                    [pth nme ext]=fileparts(fns{ind});
                    fns_dest{ind}=fullfile(pth,[sprintf('%s-%d',nme,suffix) ext]);
                    suffix=suffix+1;
                    fns_dest_full{ind}=fullfile(dest,fns_dest{ind});
                end;
                
                % Create full path
                fns_dest_full{ind}=fullfile(dest,fns_dest{ind});
            end;
            
            aas_log(aap,false,sprintf(' retrieve remote stream %s from %s:%s to %s',streamname,inputstream.host,src,dest),aap.gui_controls.colours.inputstreams);
            
            oldpth='';
            
            if (~doneremotefetch)
                % rsync in chunks. It will then compress.
                chunksize=64;
                transfernow=false;
                numtotransfer=0;
                inps='';
                for ind=1:length(fns)
                    % Copy file
                    [pth nme ext]=fileparts(fns_dest_full{ind});
                    newpth=pth;
                    if (~strcmp(oldpth,newpth))
                        aas_makedir(aap,newpth);
                        if (numtotransfer>0)
                            aas_copyfromremote(aap, inputstream.host, inps,oldpth);
                        end;
                        inps=[fullfile(src,fns{ind}) ' '];
                        oldpth=newpth;
                        numtotransfer=1;
                    else
                        inps=[inps fullfile(src,fns{ind}) ' '];
                        numtotransfer=numtotransfer+1;
                    end;
                    if (wasnamechange)
                        aas_copyfromremote(aap, inputstream.host, fullfile(src,fns{ind}),fns_dest_full{ind});
                        numtotransfer=0;
                        inps='';
                    else
                        if (numtotransfer>0) && (ind==length(fns) || numtotransfer>chunksize)
                            aas_copyfromremote(aap, inputstream.host, inps,oldpth);
                            numtotransfer=0;
                            inps='';
                        end;
                    end;
                end;
            end;
            
            % Get read to write the stream file
            [aap datecheck_md5_recalc]=aas_md5(aap,fns_dest_full,[],'filestats');
            if exist(inputstreamdesc,'file')
                delete(inputstreamdesc);
            end;
            fid_inp=fopen(inputstreamdesc,'w');
            fprintf(fid_inp,'MD5\t%s\t%s\n',md5,datecheck_md5_recalc);
            
            for ind=1:length(fns)
                % Write to stream file
                fprintf(fid_inp,'%s\n',fns_dest{ind});
            end;
            if isempty(fns)
                aas_log(aap,false,sprintf('No inputs in stream %s',streamname));
            end;
            fclose(fid_inp);
            
            gotinputs=[gotinputs;fns_dest_full];
        end;
        
        %     if (~doneremotefetch)
        %         fid=fopen(remotefetchfn,'w');
        %         fprintf(fid,sprintf('%f',now));
        %         fclose(fid);
        %     end;
    else
        %% LOCALLY - local disc or S3
        % Not remote pull
        sourcestagetag=aas_getstagetag(aap,sourcestagenumber);
        switch(aap.directory_conventions.remotefilesystem)
            case 'none'                
                % Source and destination directories
                src= aas_getpath_bydomain(aap,domain,indices,sourcestagenumber);
                dest=aas_getpath_bydomain(aap,domain,indices);
                
                outputstreamdesc=fullfile(src,sprintf('stream_%s_outputfrom_%s.txt',fromstreamname,sourcestagetag));
                
                % This is written by this routine - the input to the stage we're
                % working on
                
                % If stream will not be qualified, delete non-qualified stream
                % name, if this exists, as this will override a qualified filename, which is dangerous
                non_qualified_fn=fullfile(dest,sprintf('stream_%s_inputto_%s.txt',fromstreamname,aap.tasklist.currenttask.name));
                if ~strcmp(streamname,fromstreamname) && exist(non_qualified_fn,'file')
                    delete(non_qualified_fn);
                end;
                
                % Now produce either qualified or non-qualified (depending on
                % streamname)
                inputstreamdesc=fullfile(dest,sprintf('stream_%s_inputto_%s.txt',streamname,aap.tasklist.currenttask.name));
                
                % Copy from one drive to another on local filesystem
                %  Note the src corresponds to the outputstream of the previous
                %  stage, and the dest corresponds to the inputstream of the
                %  next stage. A touch confusing!
                if exist(outputstreamdesc,'file')
                    % Make sure output directory exists
                    [s m mid]=mkdir(dest);
                    
                    reloadfiles=true;
                    fid=fopen(outputstreamdesc,'r');
                    
                    % Load checksum with error checking
                    [aap md5]=loadmd5(aap,fid,streamname);
                    
                    % Get filenames
                    fns=textscan(fid,'%s');
                    fns=fns{1};
                    fns_dest=cell(length(fns),1);
                    fns_dest_full=cell(length(fns),1);
                    for ind=1:length(fns)
                        % Check to see whether a filename with this name has
                        % already been loaded. If so, add unique suffix
                        fns_dest{ind}=fns{ind};
                        fns_dest_full{ind}=fullfile(dest,fns_dest{ind});
                        suffix=1;
                        while (1)
                            pos=[strcmp(fns_dest_full{ind},gotinputs)];
                            if (~any(pos))
                                break;
                            end;
                            [pth nme ext]=fileparts(fns{ind});
                            fns_dest{ind}=fullfile(pth,[sprintf('%s-%d',nme,suffix) ext]);
                            suffix=suffix+1;
                            fns_dest_full{ind}=fullfile(dest,fns_dest{ind});
                        end;
                    end;
                    
                    % Also check last modified dates and file size
                    [aap datecheck_md5_recalc]=aas_md5(aap,fns_dest_full,[],'filestats');
                    % Check to see if there is already a stream file with the
                    % appropriate name in the destination, which has a matching
                    % MD5. If so, there won't be any need to copy
                    
                    if exist(inputstreamdesc,'file')
                        % Check MD5s are the same across new and previous
                        % input stream file
                        fid_out=fopen(outputstreamdesc,'r');
                        md5_lne_out=fgetl(fid_out);
                        fid_inp=fopen(inputstreamdesc,'r');
                        md5_lne=fgetl(fid_inp);
                        
                        % Check that checksums in the two files are the same
                        [junk rem]=strtok(md5_lne_out);
                        md5_o=strtok(rem);
                        [junk rem]=strtok(md5_lne);
                        md5_i=strtok(rem);
                        
                        
                        if ~strcmp(md5_o,md5_i)
                            aas_log(aap,false,sprintf('MD5 lines of stream files do not match, will recopy.'));
                        else
                            
                            % Check same number of lines in input and
                            % output files, in case copying failed half way
                            % through
                            filematch=true;
                            while ~feof(fid_out)
                                file_o=fgetl(fid_out);
                                file_i=fgetl(fid_inp);
                                if ~strcmp(deblank(file_o),deblank(file_i))
                                    filematch=false;
                                end;
                            end;
                            
                            if ~filematch
                                aas_log(aap,false,'Previous copying of files did not complete, will recopy');
                            else
                                if (length(md5_lne)>3 && (strcmp(md5_lne(1:3),'MD5')))
                                    [aap md5_inp datecheck]=loadmd5(aap,md5_lne,streamname);
                                    if (strcmp(md5_inp,md5) && strcmp(datecheck,datecheck_md5_recalc))
                                        reloadfiles=false;
                                    end;
                                    %                    fprintf('Loaded datecheck was %s and calc %s\n',datecheck,datecheck_md5_recalc);
                                end;
                            end;
                        end;
                        fclose(fid_inp);
                        fclose(fid_out);
                        
                    end;
                    
                    fclose(fid);
                    
                    if (~reloadfiles)
                        aas_log(aap,false,sprintf(' retrieve stream %s [checksum match, not recopied] from %s to %s',streamname,src,dest),aap.gui_controls.colours.inputstreams);
                    else
                        aas_log(aap,false,sprintf(' retrieve stream %s from %s to %s',streamname,src,dest),aap.gui_controls.colours.inputstreams);
                        
                        oldpth='';
                        % Get read to write the stream file
                        fid_inp=fopen(inputstreamdesc,'w');
                        fprintf(fid,'MD5\t%s\t%s\n',md5,datecheck_md5_recalc);
                        
                        
                        
                        for ind=1:length(fns)
                            % Copy file
                            [pth nme ext]=fileparts(fns_dest_full{ind});
                            newpth=pth;
                            if (~strcmp(oldpth,newpth))
                                aas_makedir(aap,newpth);
                                %relativepath_src_from_dest=relativepath(src,newpth);
                                oldpth=newpth;
                            end;
                            if (ismodified)
                                cmd=['cd ' src '; rsync -t ' fns{ind} ' ' fns_dest_full{ind}];
                            else
                                % This is a hard link, not a symlink. This
                                % takes the timestamp of the destination file,
                                % and won't be deleted if the destination is
                                % deleted. So, more like a copy...
                                cmd=['ln -f ' fullfile(src,fns{ind}) ' ' fns_dest_full{ind}];
                            end;
                            aas_shell(cmd);
                            
                            % Write to stream file
                            fprintf(fid_inp,'%s\n',fns_dest{ind});
                        end;
                        if isempty(fns)
                            aas_log(aap,false,sprintf('No inputs in stream %s',streamname));
                        end;
                        fclose(fid_inp);
                    end;
                    gotinputs=[gotinputs;fns_dest_full];
                end;
                
                
            case 's3'
                src= aas_getpath_bydomain(aap,domain,indices,'s3',sourcestagenumber);
                dest=aas_getpath_bydomain(aap,domain,indices,'s3');
                
                % Stream .txt file goes in this directory
                descriptor=sprintf('stream_%s_outputfrom_%s.txt',streamname,sourcestagetag);
                tmp=fullfile(dest,descriptor);
                
                streamentry=[aaworker.bucket '|' src '/' descriptor];
                fprintf('Trying for %s\n',streamentry);
                [aap resp]=sdb_get_attributes(aap,aaworker.streamtablename,streamentry);
                
                if (~isempty(resp))
                    fprintf('Found SDB\n');
                    s3_copyfrom_filelist(aap,dest,descriptor,aaworker.bucket,src);
                else
                    fprintf('No SDB entry found\n');
                end;
                if (exist(tmp,'file'))
                    fprintf('Found file %s\n',tmp);
                else
                    fprintf('Trying for %s\n',tmp);
                end;
                % change the name as it is now an input...
                descriptor=sprintf('stream_%s_inputto_%s.txt',streamname,aap.tasklist.currenttask.name);
                
                if (exist(tmp,'file'))
                    % prepare to fetch from AWS S3
                    [s m mid]=mkdir(dest);
                    
                    % read stream file to get input properties
                    streamdesc=fullfile(dest,descriptor);
                    movefile(tmp,streamdesc);
                    fid=fopen(streamdesc,'r');
                    
                    % Load md5 checksum from stream file with error checking
                    [aap md5stored]=loadmd5(aap,fid,streamname);
                    
                    % Get filelist from stream
                    fns=textscan(fid,'%s');
                    fns=fns{1};
                    fclose(fid);
                    
                    % Retrieve these files from the cold store
                    aas_log(aap,false,sprintf(' retrieve stream %s from s3:%s to %s',streamname,src,dest));
                    s3_copyfrom_filelist(aap,dest,fns,aaworker.bucket,src);
                    
                    % Calculate md5 of retrieved data
                    [aap md5calc]=aas_md5(aap,fns,dest);
                    
                    % Compare md5s
                    if (~strcmp(md5calc,md5stored))
                        aas_log(aap,true,sprintf('MD5 checksum does not match that of data loaded for stream %s',streamname));
                    else
                        aas_log(aap,false,' ...done, MD5 match');
                    end;
                    
                    % Check to see if any need gunzip uncompressing, if so do it
                    forgz=regexp(fns,'.*\.gz');
                    forgz=find([forgz{:}]);
                    for gzitem=forgz
                        cmd=['gunzip ' fns{gzitem}];
                        fprintf(cmd);
                        aas_shell(cmd);
                        fns{gzitem}=fns{gzitem}(1:end-3);
                    end;
                    
                    gotinputs=[gotinputs;fns];
                end;
        end;
    end;
    
end


function [aap md5 datecheck]=loadmd5(aap,fid,streamname)
if (ischar(fid))
    lne=fid;
else
    lne=fgetl(fid);
end;
pos=find(lne==9);
if (length(lne)<3 || ~strcmp(lne(1:3),'MD5') || isempty(pos))
    aas_log(aap,true,sprintf('MD5 in file %s corrupted',streamname));
end;
if (length(pos)==1)
    md5=deblank(lne(pos+1:end));
    datecheck='';
else
    md5=deblank(lne(pos(1)+1:pos(2)-1));
    datecheck=deblank(lne(pos(2)+1:end));
end;