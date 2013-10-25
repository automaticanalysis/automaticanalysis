% Make a movie of an epi time series
% [aap,resp]=aamod_movie(aap,task,i,j)
% J Carlin 20120911

function [aap,resp]=aamod_movie(aap,task,i,j)

resp='';

switch task
    case 'description'
        resp='Run aamod_movie';
    case 'summary'
        resp='Make EPI movies\n';
    case 'doit'
        for j = aap.acq_details.selected_sessions
            sesspath=aas_getsesspath(aap,i,j);
            aas_makedir(aap,sesspath);
            outfile = fullfile(sesspath,'diagnostic_aamod_movie.avi');
            % get files in this directory
            imgs=aas_getimages_bystream(aap,i,j,'epi');
            V = spm_vol(imgs);
            volumes2movie(V,outfile,aap.tasklist.currenttask.settings.fps);
        end
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end;
