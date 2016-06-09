% AA module - tsdiffana - tool to assess time series variance
% [aap,resp]=aamod_tsdiffana(aap,task,i,j)
% Rhodri Cusack MRC CBU Cambridge Aug 2004
% Tibor Auer MRC CBU Cambridge 2012-2013

function [aap,resp]=aamod_tsdiffana(aap,task,subjInd,sessInd)

resp='';

switch task
    case 'report' % Updated [TA]
        sesspath = aas_getsesspath(aap, subjInd, sessInd);
        if ~exist(fullfile(sesspath,'diagnostic_aamod_tsdiffana.jpg'),'file')
            diag(aap,subjInd,sessInd);
        end
        aap = aas_report_add(aap,subjInd,'<table><tr><td>');
        aap=aas_report_addimage(aap,subjInd,fullfile(sesspath, 'diagnostic_aamod_tsdiffana.jpg'));
        aap = aas_report_add(aap,subjInd,'</td></tr></table>');
    case 'doit'
        sesspath=aas_getsesspath(aap,subjInd,sessInd);
        
        aas_makedir(aap,sesspath);
        
        % get files in this directory
        job.imgs{1}{1} = aas_getimages_bystream(aap,subjInd,sessInd,'epi');
        job.vf = 0;
        run_tsdiffana('run','timediff',job);
        
        subjpth=aas_getsesspath(aap,subjInd,sessInd);
        aap=aas_desc_outputs(aap,subjInd,sessInd,'tsdiffana',fullfile(subjpth,'timediff.mat'));

        diag(aap,subjInd,sessInd);
        
        subjname = aas_prepare_diagnostic(aap, subjInd);
        set(gcf,'Renderer','zbuffer');
        print('-djpeg','-r150',fullfile(aap.acq_details.root, 'diagnostics', ...
            [mfilename '_' subjname '_' aap.acq_details.sessions(sessInd).name '.jpeg']));        
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end;
end

function diag(aap,subjInd,sessInd)
tsfn = aas_getfiles_bystream(aap,subjInd,sessInd,'tsdiffana');
tsdiffplot(tsfn);
try f = spm_figure('FindWin', 'Graphics'); catch; f = figure(1); end;
set(f,'Renderer','zbuffer');
print(f,'-djpeg','-r150',fullfile(aas_getsesspath(aap, subjInd, sessInd),'diagnostic_aamod_tsdiffana'));
end