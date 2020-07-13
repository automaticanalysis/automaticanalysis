% Copies across the input streams to a module from previous stages.
%
% In the case that archiving to Glacier is active, then streams that have
% been archived will have been queued for retrieval in part 1. This script
% will start working on them as they become available, and will block until
% all are retrieved and processed.
%


function [aap]=aas_retrieve_inputs_part2(aap,streamfiles)

nstreams=length(streamfiles);

% Dependencies are processed as they become available, not necessarily in
% list order
depnotdone=true(nstreams,1);

archivemessage=false(nstreams,1);

while any(depnotdone)
    for depind=1:nstreams
        if isempty(streamfiles(depind).outputstreamdesc)
            depnotdone(depind)=false;
        end;
        if depnotdone(depind)
            if (~streamfiles(depind).reloadfiles)
                aas_log(aap,false,sprintf(' retrieve stream %s [checksum match, not recopied] from %s to %s',streamfiles(depind).streamname,streamfiles(depind).src,streamfiles(depind).dest),aap.gui_controls.colours.inputstreams);
                depnotdone(depind)=false;
            else
                % Check if this stream is back from archive, if it was ever
                % there
                switch streamfiles(depind).streamlocation
                    case 'remote'
                        fid_check=fopen(streamfiles(depind).remoteoutputstreamdesc,'r');
                        line=fgetl(fid_check);
                        fclose(fid_check);
                    case 'local'
                        fid_check=fopen(streamfiles(depind).outputstreamdesc,'r');
                        line=fgetl(fid_check);
                        fclose(fid_check);
                end;
                
                % If not archived
                if strcmp(line(1:11),'ARCHIVED TO')
                    if ~archivemessage(depind)
                        aas_log(aap,false,sprintf(' stream %s is archived, awaiting retrieval',streamfiles(depind).streamname),aap.gui_controls.colours.inputstreams);
                        archivemessage(depind)=true;
                    end;
                else
                    % Retrieve remote or local stream?
                    switch streamfiles(depind).streamlocation
                        case 'remote'
                            aas_log(aap,false,sprintf(' retrieve remote stream %s from %s:%s to %s',streamfiles(depind).streamname,streamfiles(depind).inputstream.host,streamfiles(depind).src,streamfiles(depind).dest),aap.gui_controls.colours.inputstreams);
                            
                            oldpth='';
                            
                            % rsync in chunks. It will then compress.
                            chunksize=64;
                            transfernow=false;
                            numtotransfer=0;
                            inps='';
                            for ind=1:length(streamfiles(depind).fns)
                                % Copy file
                                [pth,nme,ext]=fileparts(streamfiles(depind).fns_dest_full{ind});
                                newpth=pth;
                                if (~strcmp(oldpth,newpth))
                                    aas_makedir(aap,newpth);
                                    if (numtotransfer>0)
                                        aas_copyfromremote(aap, streamfiles(depind).inputstream.host, inps,oldpth,'verbose',0,'allowcache',streamfiles(depind).inputstream.allowcache);
                                    end;
                                    inps=[fullfile(streamfiles(depind).src,streamfiles(depind).fns{ind}) ' '];
                                    oldpth=newpth;
                                    numtotransfer=1;
                                else
                                    inps=[inps fullfile(streamfiles(depind).src,streamfiles(depind).fns{ind}) ' '];
                                    numtotransfer=numtotransfer+1;
                                end;
                                if (streamfiles(depind).wasnamechange)
                                    aas_copyfromremote(aap, streamfiles(depind).inputstream.host, fullfile(streamfiles(depind).src,streamfiles(depind).fns{ind}),streamfiles(depind).fns_dest_full{ind},'verbose',0,'allowcache',streamfiles(depind).inputstream.allowcache);
                                    numtotransfer=0;
                                    inps='';
                                else
                                    if (numtotransfer>0) && (ind==length(streamfiles(depind).fns) || numtotransfer>chunksize)
                                        aas_copyfromremote(aap, streamfiles(depind).inputstream.host, inps,oldpth,'verbose',0,'allowcache',streamfiles(depind).inputstream.allowcache);
                                        numtotransfer=0;
                                        inps='';
                                    end;
                                end;
                            end;
                            
                            % Get read to write the stream file
                            [aap,datecheck_md5_recalc]=aas_md5(aap,streamfiles(depind).fns_dest_full,[],'filestats');
                            if exist(streamfiles(depind).inputstreamdesc,'file')
                                try
                                    delete(streamfiles(depind).inputstreamdesc);
                                catch
                                end;
                            end;
                            fid_inp=fopen(streamfiles(depind).inputstreamdesc,'w');
                            fprintf(fid_inp,'MD5\t%s\t%s\n',streamfiles(depind).md5,datecheck_md5_recalc);
                            
                            for ind=1:length(streamfiles(depind).fns)
                                % Write to stream file
                                fprintf(fid_inp,'%s\n',streamfiles(depind).fns_dest{ind});
                            end;
                            if isempty(streamfiles(depind).fns)
                                aas_log(aap,false,sprintf('No inputs in stream %s',streamname));
                            end;
                            fclose(fid_inp);
                            depnotdone(depind)=false;
                            
                        case 'local'
                            
                            aas_log(aap,false,sprintf(' retrieving stream %s from %s to %s',streamfiles(depind).streamname,streamfiles(depind).src,streamfiles(depind).dest),aap.gui_controls.colours.inputstreams);
                            oldpth='';
                            % Get ready to write the stream file
                            fid_inp=fopen(streamfiles(depind).inputstreamdesc,'w');
                            fprintf(fid_inp,'MD5\t%s\t%s\n',streamfiles(depind).md5,streamfiles(depind).datecheck_md5_recalc);
                            for ind=1:length(streamfiles(depind).fns_dest_full)
                                % Copy file
                                [pth nme ext]=fileparts(streamfiles(depind).fns_dest_full{ind});
                                newpth=pth;
                                if (~strcmp(oldpth,newpth))
                                    aas_makedir(aap,newpth);
                                    %relativepath_src_from_dest=relativepath(src,newpth);
                                    oldpth=newpth;
                                end;
                                if (streamfiles(depind).ismodified) || ~aap.options.hardlinks
                                    cmd=['cd ' streamfiles(depind).src '; rsync -t ' streamfiles(depind).fns{ind} ' ' streamfiles(depind).fns_dest_full{ind}];
                                else
                                    % This is a hard link, not a symlink. This
                                    % takes the timestamp of the destination file,
                                    % and won't be deleted if the destination is
                                    % deleted. So, more like a copy...
                                    cmd=['ln -f ' fullfile(streamfiles(depind).src,streamfiles(depind).fns{ind}) ' ' streamfiles(depind).fns_dest_full{ind}];
                                end;
                                aas_shell(cmd);
                                
                                % Write to stream file
                                fprintf(fid_inp,'%s\n',streamfiles(depind).fns_dest{ind});
                            end;
                            
                            
                            % Alert archiving system to new data if it is at work here
                            [spth snme sext]=fileparts(streamfiles(depind).inputstreamdesc);
                            dotpath=fullfile(streamfiles(depind).dest,['.' snme]);
                            
                            if exist(dotpath,'dir')
                                localchangelog=fullfile(dotpath,'log_localchanges.txt');
                                lclfid=fopen(localchangelog,'a');
                                fprintf(lclfid,'%s\tSTREAM REWRITTEN BY aa\n',datestr(now,'yyyy-mm-ddTHH:MM:SS.FFF'));
                                fclose(lclfid);
                            end;
                            
                            %                             % Bring across archiving if present, adjusting
                            %                             % as necessary to ensure duplicate (hard
                            %                             % linked) streams are archived only once, but
                            %                             % deleted twice
                            %                             [pth nme ext]=fileparts(streamfiles(depind).outputstreamdesc);
                            %                             srcdotpath=fullfile(streamfiles(depind).src,['.' nme],[streamfiles(depind).md5 'v1']);
                            %                             [pth nme ext]=fileparts(streamfiles(depind).inputstreamdesc);
                            %                             destdotpath=fullfile(streamfiles(depind).dest,['.' nme],[streamfiles(depind).md5 'v1']);
                            %                             if exist(srcdotpath,'dir')
                            %                                 aas_log(aap,false,sprintf('Found archive %s',srcdotpath));
                            %                                 if ~exist(destdotpath,'dir')
                            %                                     aas_makedir(aap,destdotpath);
                            %                                     % Always just copy local change log, as each end
                            %                                     % has its own independent life
                            %                                     srclocalchangelog=fullfile(srcdotpath,'log_localchanges.txt');
                            %                                     if exist(srclocalchangelog,'file')
                            %                                         copyfile(srclocalchangelog,destdotpath);
                            %                                     end;
                            %                                     % But symlink s3 log if not modified, so we don't
                            %                                     % duplicate backing up
                            %                                     srcs3changelog=fullfile(srcdotpath,'log_s3changes.txt');
                            %                                     if (streamfiles(depind).ismodified)
                            %                                         copyfile(srcs3changelog,destdotpath);
                            %                                     else
                            %                                         aas_shell(['ln -f ' srcs3changelog ' ' destdotpath]);
                            %                                     end;
                            %                                 end;
                            %                             end;
                            
                            if isempty(streamfiles(depind).fns)
                                aas_log(aap,false,sprintf('No inputs in stream %s',streamfiles(depind).streamname));
                            end;
                            fclose(fid_inp);
                            
                            depnotdone(depind)=false;
                    end;
                end;
            end;
        end;
    end;
    
    % Wait a moment to avoid overloading the filesystem
    pause(1.0);
end