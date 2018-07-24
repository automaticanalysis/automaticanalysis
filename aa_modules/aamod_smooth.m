% AA module - smoothing
% [aap,resp]=aamod_smooth(aap,task,subj,sess)
% Gaussian smoothing images, important for Gaussian Field Theory
% Kernel size determined by aap.spm_analysis.FWHM
% Rhodri Cusack MRC CBU Cambridge Jan 2006- Aug 2007
% Now once per session for improved performance when running in parallel

function [aap,resp]=aamod_smooth(aap,task,subj,sess)
resp='';

switch task
    case 'report'
        
    case 'doit'
        
        % Is session specified in task header?
        if (isfield(aap.tasklist.currenttask.settings,'session'))
            sess = aap.tasklist.currenttask.settings.session;
        end
        
        streams=aap.tasklist.currenttask.inputstreams.stream;
        
        for streamind=1:length(streams)
            
            % Images to smooth
            if (exist('sess','var'))
                P = aas_getfiles_bystream(aap,aap.tasklist.currenttask.domain,[subj sess],streams{streamind});
            else
                P = aas_getfiles_bystream(aap,subj,streams{streamind});
            end
            
            % now smooth
            s   = aap.tasklist.currenttask.settings.FWHM;
            n   = size(P,1);
            
            spm_progress_bar('Init',n,'Smoothing','Volumes Complete');
            outputfns=[];
            for imnum = 1:n
                Q = deblank(P(imnum,:));
                [pth,nam,xt,nm] = spm_fileparts(deblank(Q));
                fn=['s' nam xt nm];
                U = fullfile(pth,fn);
                outputfns=strvcat(outputfns,U);
                % Ignore .hdr files from this list...
                if isempty(strfind(P(n,:), '.hdr'))
                    spm_smooth(Q,U,s);
                    spm_progress_bar('Set',imnum);
                end
            end
            
            % Describe outputs
            if (exist('sess','var'))
                aap=aas_desc_outputs(aap,aap.tasklist.currenttask.domain,[subj sess],streams{streamind},outputfns);
            else
                aap=aas_desc_outputs(aap,subj,streams{streamind},outputfns);
            end
            
        end;
        % All done
        spm_progress_bar('Clear');
    case 'checkrequirements'
        in =  aas_getstreams(aap,'input');
        [stagename, index] = strtok_ptrn(aap.tasklist.currenttask.name,'_0');
        stageindex = sscanf(index,'_%05d');
        out = aap.tasksettings.(stagename)(stageindex).outputstreams.stream; if ~iscell(out), out = {out}; end
        for s = 1:numel(in)
            if ~strcmp(out{s},in{s})
                aap = aas_renamestream(aap,aap.tasklist.currenttask.name,out{s},in{s},'output');
            end
            aas_log(aap,false,['INFO: ' aap.tasklist.currenttask.name ' output stream: ''' in{s} '''']);
        end        
end



