% AA module - tsdiffana - tool to assess time series variance
% [aap,resp]=aamod_tsdiffana(aap,task,i,j)
% Rhodri Cusack MRC CBU Cambridge Aug 2004

function [aap,resp]=aamod_tsdiffana(aap,task,i,j)

resp='';

switch task
    case 'domain'
        resp='session';   % this module needs to be run once per session
    case 'whentorun'
        resp='justonce';  % should this be run everytime or justonce?        
    case 'description'
        resp='Run tsdiffana';
    case 'summary'
        resp='Check time series variance using tsdiffana\n';
    case 'report'
        aap.report.html=strcat(aap.report.html,'<table><tr><td>');
        aap=aas_report_addimage(aap,fullfile(aas_getsesspath(aap,i,j),'diagnostic_aamod_tsdiffana.jpg'));     
        aap.report.html=strcat(aap.report.html,'</td></tr></table>');
    case 'doit'
        subjpath=aas_getsubjpath(aap,i);
        sesspath=aas_getsesspath(aap,i,j);

        aas_makedir(aap,sesspath);
        
        % imgs=spm_get('files',sesspath,[aap.directory_conventions.subject_filenames{i} '*nii']); % changed img to nii [djm160206]
 
        % added in place of previous line [djm 160206]...
            % get the subdirectories in the main directory
            dirn = aas_getsesspath(aap,i,j);
            % get files in this directory
            imgs=aas_getimages_bystream(aap,i,j,'epi');
            
        tsdiffana(imgs,0);
        
        % Now produce graphical check
        figure(spm_figure('FindWin', 'Graphics'));
        print('-djpeg','-r75',fullfile(sesspath,'diagnostic_aamod_tsdiffana'));
     
    
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end;