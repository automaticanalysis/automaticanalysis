% AA module - IMPORT FILES AS STREAM from NIFTI

function [aap,resp]=aamod_importfilesasstream(aap,task,varargin)

resp='';

switch task
    case 'report'
    case 'doit'
        streamname=aap.tasklist.currenttask.outputstreams.stream{1};
        
        % Find where this module is working
        switch(length(varargin))
            case 0
                currpth=aas_getstudypath(aap);
            case 1
                currpth=aas_getsubjpath(aap,varargin{1});
            case 2
                currpth=aas_getsesspath(aap,varargin{1},varargin{2});
        end;

        % Copy the files to import into the current module's directory
        allnames={};
        localnames={};
        for fileind=1:length(aap.tasklist.currenttask.settings.filestoimport)
            [pth nme ext]=fileparts(aap.tasklist.currenttask.settings.filestoimport{fileind});
            timesalready=sum([strcmp([nme ext],allnames)]);
            allnames{end+1}=[nme ext];
            localnames{fileind}=fullfile(currpth,sprintf('%s-%d.%s',nme,timesalready,ext));
            copyfile(aap.tasklist.currenttask.settings.filestoimport{fileind},localnames{fileind});
        end;
    
        % Describe the output stream
        switch(length(varargin))
            case 0
                aap=aas_desc_outputs(aap,streamname,localnames);
            case 1
                aap=aas_desc_outputs(aap,varargin{1},streamname,localnames);
            case 2
                aap=aas_desc_outputs(aap,varargin{1},varargin{2},streamname,localnames);
        end;

    case 'checkrequirements'
        streamname=aap.tasklist.currenttask.outputstreams.stream{1};
        fti=aap.tasklist.currenttask.settings.filestoimport;
        if (iscell(fti))
            for fileind=1:length(fti)
                if (~exist(fti{fileind},'file'))
                    aas_log(aap,true,sprintf('File marked as initial import for stream %s not found: %s',streamname,fti{fileind}));
                end;
            end;
        else
            for fileind=1:size(fti,1)
                if (~exist(fti(fileind,:),'file'))
                    aas_log(aap,true,sprintf('File marked as initial import for stream %s not found: %s',streamname,fti(fileind,:)));
                end;
            end;
        end;
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end;














