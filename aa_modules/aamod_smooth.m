% AA module - smoothing
% [aap,resp]=aamod_smooth(aap,task,i,j)
% Gaussian smoothing images, important for Gaussian Field Theory
% Kernel size determined by aap.spm_analysis.FWHM
% Rhodri Cusack MRC CBU Cambridge Jan 2006- Aug 2007
% Now once per session for improved performance when running in parallel

function [aap,resp]=aamod_smooth(aap,task,i,j)
resp='';

switch task
    case 'domain'
        resp='session';   % this module needs to be run once per session
        
    case 'description'
        resp='SPM5 smooth';
        
    case 'summary'
        subjpath=aas_getsubjpath(i);
        resp=sprintf('Smooth run on %s\n',subjpath);
    case 'report'
        
    case 'doit'
                
        streams=aap.tasklist.currenttask.inputstreams.stream;
        
        for streamind=1:length(streams)
            
            subj.imgs=aas_getimages_bystream(aap,i,j,streams{streamind});
            
            % now smooth
            s   = aap.tasklist.currenttask.settings.FWHM;
            n   = size(subj.imgs,1);
            
            spm_progress_bar('Init',n,'Smoothing','Volumes Complete');
            outputfns=[];
            for imnum = 1:n
                Q = deblank(subj.imgs(imnum,:));
                [pth,nam,xt,nm] = spm_fileparts(deblank(Q));
                fn=['s' nam xt nm];
                U = fullfile(pth,fn);
                outputfns=strvcat(outputfns,U);
                spm_smooth(Q,U,s);
                spm_progress_bar('Set',imnum);
            end
            
            % Describe outputs
            aap=aas_desc_outputs(aap,i,j,streams{streamind},outputfns);
            
        end;
        % All done
        spm_progress_bar('Clear');
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end;



