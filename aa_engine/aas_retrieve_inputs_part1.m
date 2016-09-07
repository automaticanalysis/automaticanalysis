% Retrieves inputs from a particular stream to a particular destination
%

function [gotinputs streamfiles]=aas_retrieve_inputs_part1(aap,inputstream,gotinputs,deps)
global aaworker

streamfiles=[];

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
        [junk,isCreated] = aas_makedir(aap,dest);        
        
        % First we need to get the previous aap file
        remote_aap_fn=fullfile(dest,'remote_aap_parameters.mat');
        
        aap=aas_copyfromremote(aap,inputstream.host,inputstream.aapfilename,remote_aap_fn,'allowcache',inputstream.allowcache,'verbose', 0);
        
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
        aap_remote=aas_setcurrenttask(aap_remote,modind,'nodefault');
        
        % Check for recent additions to aap in case we're using an older
        % version on the remote machine
        if ~isfield(aap_remote.directory_conventions,'parallel_dependencies')
            aap_remote.directory_conventions.parallel_dependencies=aap.directory_conventions.parallel_dependencies;
            aas_log(aap,false,'Warning: Remote aap does not contain parallel_dependencies field - probably from older aa');
        end;
        
        % Map the local indices to remote indices
        indexMap = aas_mapindices_betweenAAPs(aap, aap_remote);
        localTree = aas_dependencytree_finddomain(domain, aap.directory_conventions.parallel_dependencies, {});
        remotespecified = arrayfun(@(x,y) ~isempty(indexMap.(x{1})), localTree); % backward compatibility for pipelines with earlier aa version
        remoteindices = arrayfun(@(x,y) indexMap.(x{1})(y), localTree(remotespecified), [1 indices(remotespecified(2:end))]);
        remoteindices(1) = []; % pop the study index, we never seem to use that.
        remoteindices(~remotespecified(2:end)) = 0;
        
        if any(remoteindices == 0)
            badIndex = find(remoteindices==0,1,'first');
            domainItems = aas_getNames_bydomain(aap, localTree{badIndex+1});
            aas_log(aap, false, sprintf('WARNING: Remote AAP doesn''t have %s ''%s''!', localTree{badIndex+1}, domainItems{1}{indices(badIndex)}));
        end
                
        % IT's possible that AAPs from different versions of AA have
        % different dependency lists... let's not error if that;s the case.
        %  Instead, just don't bother searching for it.
        if ~isempty(aas_dependencytree_finddomain(domain, aap_remote.directory_conventions.parallel_dependencies,{}))
            
            try src=aas_getpath_bydomain(aap_remote,domain,remoteindices,modind);
            catch, continue; end % skip domain if not exist in remote
            
            remoteoutputstreamdesc=fullfile(src,sprintf('stream_%s_outputfrom_%s.txt',fromstreamname,inputstream.sourcestagename));
            
            % This is written by this routine - the input to the stage we're
            % working on
            
            % Delete non-qualified stream name, if this exists, as this will
            % override a qualified filename, which is dangerous
            %             non_qualified_fn=fullfile(dest,sprintf('stream_%s_inputto_%s.txt',fromstreamname,aap.tasklist.currenttask.name));
            %             if (exist(non_qualified_fn,'file'))
            %                 delete(non_qualified_fn);
            %             end;
            
            % Now produce either qualified or non-qualified (depending on
            % streamname)
            inputstreamdesc=fullfile(dest,sprintf('stream_%s_inputto_%s.txt',streamname,aap.tasklist.currenttask.name));
            
            % Copy from one drive to another on local filesystem
            %  Note the src corresponds to the outputstream of the previous
            %  stage, and the dest corresponds to the inputstream of the
            %  next stage. A touch confusing!
            outputstreamdesc=fullfile(dest,sprintf('stream_%s_remoteoutputfrom_%s_%s.txt',fromstreamname,inputstream.host,inputstream.sourcestagename));
            
            aap=aas_copyfromremote(aap,inputstream.host,remoteoutputstreamdesc,outputstreamdesc,'allow404',1,'allowcache',inputstream.allowcache, 'verbose', 1);
            
            if (exist(outputstreamdesc,'file'))
                
                % Register access to this stream file
                aas_archive_register_access(aap,outputstreamdesc);
                
                reloadfiles=true;
                fid=fopen(outputstreamdesc,'r');
                
                % Load checksum with error checking
                [aap,md5,datecheck, isarchived]=aas_load_md5(aap,fid,streamname);

                % Issue retrieval request if this stream is archived
                if isarchived
                    % Have to retrieve at the remote location, as this is
                    % where the archiving info is stored
                    aas_log(aap,false,sprintf(' requesting retrieval of %s',streamname), aap.gui_controls.colours.inputstreams);
                    aas_archive_request_retrieval(aap,remoteoutputstreamdesc);
                end;

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
                            if ~strcmp(deblank(file_o),deblank(file_i)) || ~exist(fullfile(dest,deblank(file_i)),'file')
                                filematch=false;
                            end;
                        end;
                        
                        if ~filematch
                            aas_log(aap,false,'Previous copying of files did not complete, will recopy');
                        else
                            if (length(md5_lne)>3 && (strcmp(md5_lne(1:3),'MD5')))
                                [aap md5_inp datecheck]=aas_load_md5(aap,md5_lne,streamname);
                                if (strcmp(md5_inp,md5) && strcmp(datecheck,datecheck_md5_recalc))
                                    reloadfiles=false;
                                end;
                                aas_log(aap,false,sprintf('Loaded datecheck was %s and calc %s',datecheck,datecheck_md5_recalc));
                            end;
                        end;
                    end;
                    fclose(fid_inp);
                    fclose(fid_out);
                    
                end;
                
                fclose(fid);
                
                if ~reloadfiles
                    aas_log(aap,false,sprintf(' retrieve stream %s [checksum match, not recopied] from %s to %s',streamname,src,dest),aap.gui_controls.colours.inputstreams);
                else
                
                    % Just queue for retrieval in part 2 (remote)                   
                    streamfiles(depind).streamname=streamname;
                    streamfiles(depind).src=src;
                    streamfiles(depind).dest=dest;
                    streamfiles(depind).inputstreamdesc=inputstreamdesc;
                    streamfiles(depind).outputstreamdesc=outputstreamdesc;
                    streamfiles(depind).remoteoutputstreamdesc=remoteoutputstreamdesc;
                    streamfiles(depind).md5=md5;
                    streamfiles(depind).fns_dest_full=fns_dest_full;
                    streamfiles(depind).fns_dest=fns_dest;
                    streamfiles(depind).fns=fns;
                    streamfiles(depind).reloadfiles=reloadfiles;
                    streamfiles(depind).datecheck_md5_recalc=datecheck_md5_recalc;
                    streamfiles(depind).ismodified=ismodified;
                    streamfiles(depind).streamlocation='remote';
                    streamfiles(depind).inputstream=inputstream;
                    streamfiles(depind).wasnamechange=wasnamechange;
                end
                
                gotinputs=[gotinputs;fns_dest_full];
            else % cleanup
                if isCreated, rmdir(dest); end
            end;
            
        end;
        
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
                    % Register access to this stream file
                    aas_archive_register_access(aap,outputstreamdesc);
                    
                    % Make sure output directory exists
                    aas_makedir(aap,dest);
                    
                    reloadfiles=true;
                    fid=fopen(outputstreamdesc,'r');
                    
                    % Load checksum with error checking
                    [aap,md5,datecheck, isarchived]=aas_load_md5(aap,fid,streamname);
                    
                    % Issue retrieval request if this stream is archived
                    if isarchived
                        aas_log(aap,false,sprintf(' requesting retrieval of %s',streamname), aap.gui_controls.colours.inputstreams);
                        aas_archive_request_retrieval(aap,outputstreamdesc);
                    end;
                    
                    % Get data filenames
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
                                if feof(fid_inp)
                                    filematch=false;
                                    break;
                                end;
                                file_o=fgetl(fid_out);
                                file_i=fgetl(fid_inp);
                                if ~strcmp(deblank(file_o),deblank(file_i)) || ~exist(fullfile(dest,deblank(file_i)),'file')
                                    filematch=false;
                                end;
                            end;
                            
                            if ~filematch
                                aas_log(aap,false,'Previous copying of files did not complete, will recopy');
                            else
                                if (length(md5_lne)>3 && (strcmp(md5_lne(1:3),'MD5')))
                                    [aap md5_inp datecheck]=aas_load_md5(aap,md5_lne,streamname);
                                    if (strcmp(md5_inp,md5) && strcmp(datecheck,datecheck_md5_recalc))
                                        reloadfiles=false;
                                    end;
                                    aas_log(aap,false,sprintf('Loaded datecheck was %s and calc %s',datecheck,datecheck_md5_recalc));
                                end;
                            end;
                        end;
                        fclose(fid_inp);
                        fclose(fid_out);
                        
                    end;
                    
                    fclose(fid);
                    
                    % These are stored for part 2
                    streamfiles(depind).streamname=streamname;
                    streamfiles(depind).src=src;
                    streamfiles(depind).dest=dest;
                    streamfiles(depind).inputstreamdesc=inputstreamdesc;
                    streamfiles(depind).outputstreamdesc=outputstreamdesc;
                    streamfiles(depind).md5=md5;
                    streamfiles(depind).fns_dest_full=fns_dest_full;
                    streamfiles(depind).fns_dest=fns_dest;
                    streamfiles(depind).fns=fns;
                    streamfiles(depind).reloadfiles=reloadfiles;
                    streamfiles(depind).datecheck_md5_recalc=datecheck_md5_recalc;
                    streamfiles(depind).ismodified=ismodified;
                    streamfiles(depind).streamlocation='local';
                    
                    % Actual copying of data files now in part 2
                    % This allows asynchronous retrieval of all of the
                    % stream files together, when they've gone off to S3
                    % or Glacier
                    
                    gotinputs=[gotinputs;streamfiles(depind).fns_dest_full];
                end;
                
                
            case 's3'
                src= aas_getpath_bydomain(aap,domain,indices,'s3',sourcestagenumber);
                dest=aas_getpath_bydomain(aap,domain,indices,'s3');
                
                % Stream .txt file goes in this directory
                descriptor=sprintf('stream_%s_outputfrom_%s.txt',streamname,sourcestagetag);
                tmp=fullfile(dest,descriptor);
                
                streamentry=[aaworker.bucket '|' src '/' descriptor];
                aas_log(aap,false,sprintf('Trying for %s',streamentry));
                [aap resp]=sdb_get_attributes(aap,aaworker.streamtablename,streamentry);
                
                if (~isempty(resp))
                    aas_log(aap,false,'Found SDB');
                    s3_copyfrom_filelist(aap,dest,descriptor,aaworker.bucket,src);
                else
                    aas_log(aap,false,'No SDB entry found');
                end;
                if (exist(tmp,'file'))
                    aas_log(aap,false,sprintf('Found file %s',tmp));
                else
                    aas_log(aap,false,sprintf('Trying for %s',tmp));
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
                    [aap md5stored]=aas_load_md5(aap,fid,streamname);
                    
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
                        aas_shell(cmd);
                        fns{gzitem}=fns{gzitem}(1:end-3);
                    end;
                    
                    gotinputs=[gotinputs;fns];
                end;
        end;
    end;
    
end


