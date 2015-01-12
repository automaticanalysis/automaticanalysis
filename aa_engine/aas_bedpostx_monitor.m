function [aap alldone]=aas_bedpostx_monitor(aap,diff_pth,bedpost_pth, waitforcompletion, isfirsttime)
fslext=aas_getfslext(aap);
alldone=false;

if isfirsttime
    % Wait for diff_slices directory to appear
    while ~exist(fullfile(bedpost_pth,'diff_slices'),'dir')
        pause(5.0);
    end;
end;

% Now wait for it to disappear, while counting
oldcompleted=-1;
while(true)
    if ~exist(fullfile(bedpost_pth,'diff_slices'),'dir') && exist(fullfile(bedpost_pth,['dyads1' fslext]),'file')
        alldone=true;
        aas_log(aap,false,'bedpostX has finished!');
        break;
    else
        % Get number of slices
        cmd=['fslinfo ' diff_pth];
        [s w]=aas_runfslcommand(aap,cmd);
        completedslices=0;
        
        if ~s
            pos_start=strfind(w,'dim3');
            if ~isempty(pos_start)
                nslices=str2num(strtok(deblank(w(pos_start(1)+5:end)),10));
                
                % Now check which have finished
                
                for sliceind=0:(nslices-1)
                    if exist(fullfile(bedpost_pth,'diff_slices',sprintf('data_slice_%04d',sliceind),['dyads1' fslext]),'file')
                        completedslices=completedslices+1;
                    end;
                end;
            end;
        end;
        
        if oldcompleted~=completedslices & exist('nslices','var')
            aas_log(aap,false,sprintf(' at %s completed %d/%d for %s',datestr(now,31),completedslices,nslices,bedpost_pth));
            oldcompleted=completedslices;
        end;
    end;
    pause(10.0);
    if ~waitforcompletion
        break
    end;
end;
