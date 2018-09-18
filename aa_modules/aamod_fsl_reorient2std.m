function [aap,resp]=aamod_fsl_reorient2std(aap,task,varargin)
resp='';

switch task
    case 'report'
    case 'doit'
        % obtain input filename(s)
        indices = cell2mat(varargin);        
        streams=aas_getstreams(aap,'input');
        fname = aas_getfiles_bystream(aap,aap.tasklist.currenttask.domain,indices,streams{1});
        if size(fname,1)>1
            aas_log(aap,true,'no support for multiple inputs at present!');
        end

        % figure out output name
        [datadir,fn,ext] = fileparts(fname);
        fslext=aas_getfslext(aap);
        outname = fullfile(datadir,['reo_',fn,fslext]);

        % run it
        [junk, w]=aas_runfslcommand(aap,...
            sprintf('fslreorient2std %s %s',fname,outname));
        
        % Describe outputs
        aap=aas_desc_outputs(aap,aap.tasklist.currenttask.domain,indices,streams{1},outname);
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end
