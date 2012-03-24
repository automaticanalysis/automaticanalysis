% Retrieves inputs from a particular stream to a particular destination
%

function [gotinputs]=aas_retrieve_inputs(aap,streamname,sourcestagenumber,gotinputs,varargin)
global aaworker

% Only take part of stream after last period
pos=find(streamname=='.');
if (~isempty(pos))
    fromstreamname=streamname(pos(end)+1:end);
else
    fromstreamname=streamname;
end;

sourcestagetag=aas_getstagetag(aap,sourcestagenumber);
switch(aap.directory_conventions.remotefilesystem)
    case 'none'
        switch length(varargin)
            case 0
                src=aas_getstudypath(aap,sourcestagenumber);
                dest=aas_getstudypath(aap);
            case 1
                src=aas_getsubjpath(aap,varargin{1},sourcestagenumber);
                dest=aas_getsubjpath(aap,varargin{1});
            case 2
                src=aas_getsesspath(aap,varargin{1},varargin{2},sourcestagenumber);
                dest=aas_getsesspath(aap,varargin{1},varargin{2});
        end;
        
        outputstreamdesc=fullfile(src,sprintf('stream_%s_outputfrom_%s.txt',fromstreamname,sourcestagetag));
        
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
        if (exist(outputstreamdesc,'file'))
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
                suffix=1;
                while (1)
                    pos=[strcmp(fns_dest{ind},gotinputs)];
                    if (~any(pos))
                        break;
                    end;
                    [pth nme ext]=fileparts(fns{ind})
                    fns_dest{ind}=fullfile(pth,[sprintf('%s-%d',nme,suffix) ext]);
                    suffix=suffix+1;
                end;
                
                % Create full path
                fns_dest_full{ind}=fullfile(dest,fns_dest{ind});
            end;
            
            % Also check last modified dates and file size
            [aap datecheck_md5_recalc]=aas_md5(aap,fns_dest_full,[],'filestats');
            % Check to see if there is already a stream file with the
            % appropriate name in the destination, which has a matching
            % MD5. If so, there won't be any need to copy
            if (exist(inputstreamdesc,'file'))
                fid_inp=fopen(inputstreamdesc,'r');
                
                md5_lne=fgetl(fid_inp);
                     

                if (length(md5_lne)>3 && (strcmp(md5_lne(1:3),'MD5')))
                    [aap md5_inp datecheck]=loadmd5(aap,md5_lne,streamname);
                    if (strcmp(md5_inp,md5) && strcmp(datecheck,datecheck_md5_recalc))
                        reloadfiles=false;
                    end;
%                    fprintf('Loaded datecheck was %s and calc %s\n',datecheck,datecheck_md5_recalc);
                end;
                fclose(fid_inp);
                
                
            end;
            
            if (~reloadfiles)
                aas_log(aap,false,sprintf(' retrieve stream %s [checksum match, not recopied] from %s to %s',streamname,src,dest),aap.gui_controls.colours.inputstreams);
            else
                aas_log(aap,false,sprintf(' retrieve stream %s from %s to %s',streamname,src,dest),aap.gui_controls.colours.inputstreams);
                
                fclose(fid);
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
                        oldpth=newpth;
                    end;
                    cmd=['cd ' src '; rsync -t ' fns{ind} ' ' fns_dest_full{ind}];
                    aas_shell(cmd);
                    
                    % Write to stream file
                    fprintf(fid_inp,'%s\n',fns_dest{ind});
                end;
                if isempty(fns)
                    aas_log(aap,false,sprintf('No inputs in stream %s',streamname));
                end;
                fclose(fid_inp);
            end;
            gotinputs=fns;
        else
            gotinputs={};
        end;
        
    case 's3'
        switch length(varargin)
            case 0
                src=aas_getstudypath(aap,'s3',sourcestagenumber);
                dest=aas_getstudypath(aap);
            case 1
                src=aas_getsubjpath(aap,varargin{1},'s3',sourcestagenumber);
                dest=aas_getsubjpath(aap,varargin{1});
            case 2
                src=aas_getsesspath(aap,varargin{1},varargin{2},'s3',sourcestagenumber);
                dest=aas_getsesspath(aap,varargin{1},varargin{2});
        end;

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
            
            gotinputs=fns;
        else
            gotinputs={};
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