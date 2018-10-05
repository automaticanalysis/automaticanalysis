% AA module - tsdiffana - tool to assess time series variance
% [aap,resp]=aamod_tsdiffana(aap,task,i,j)
% Rhodri Cusack MRC CBU Cambridge Aug 2004
% Tibor Auer MRC CBU Cambridge 2012-2013

function [aap,resp]=aamod_tsdiffana(aap,task,subjInd,sessInd)

resp='';

switch task
	
    case 'report' % Updated [TA]
		
        sesspath = aas_getsesspath(aap, subjInd, sessInd);
        imgList = spm_select('FPList',sesspath,'^diagnostic_aamod_tsdiffana.*jpg$');
        if isempty(imgList)
            diag(aap,subjInd,sessInd);
            imgList = spm_select('FPList',sesspath,'^diagnostic_aamod_tsdiffana.*jpg$');
        end
        aap = aas_report_add(aap,subjInd,'<table><tr><td>');
        for n = 1:size(imgList,1)
            aap=aas_report_addimage(aap,subjInd,imgList(n,:));
        end
        aap = aas_report_add(aap,subjInd,'</td></tr></table>');
		
    case 'doit'
		
        sesspath=aas_getsesspath(aap,subjInd,sessInd);
        
        aas_makedir(aap,sesspath);
		
        job.imgs{1}{1} = aas_getfiles_bystream(aap,subjInd,sessInd,'epi');
        job.vf = 0;
        run_tsdiffana('run','timediff',job);
        
        subjpth=aas_getsesspath(aap,subjInd,sessInd);
        aap=aas_desc_outputs(aap,subjInd,sessInd,'tsdiffana',fullfile(subjpth,'timediff.mat'));

        diag(aap,subjInd,sessInd,true);
        
        subjname = aas_prepare_diagnostic(aap, subjInd);
        for i = 1:2
            f = spm_figure('FindWin', sprintf('Graphics%d',i));
            set(f,'Renderer','zbuffer');
            print(f, '-djpeg','-r150',fullfile(aap.acq_details.root, 'diagnostics', ...
                [mfilename '_' subjname '_' aap.acq_details.sessions(sessInd).name sprintf('_%d.jpeg',i)]));
            spm_figure('Close',f)
        end

    case 'checkrequirements'
        
	otherwise
		
        aas_log(aap,1,sprintf('Unknown task %s',task));
end;
end

function diag(aap,subjInd,sessInd,keepOpen)
tdfn = aas_getfiles_bystream(aap,subjInd,sessInd,'tsdiffana');

opts = {};
if aas_stream_has_contents(aap,[subjInd,sessInd],'realignment_parameter')
    mp = cellstr(aas_getfiles_bystream(aap,subjInd,sessInd,'realignment_parameter'));
    opts = {'r' mp{cell_index(mp,'.txt')}};
end

f(1) = spm_figure('Create', 'Graphics1'); spm_figure('Clear',f(1),'Graphics1');
f(2) = spm_figure('Create', 'Graphics2'); spm_figure('Clear',f(2),'Graphics2');
tsdiffplot(tdfn,f,opts{:});
for i = 1:2
    set(f(i),'Renderer','zbuffer');
    print(f(i),'-djpeg','-r150',fullfile(aas_getsesspath(aap, subjInd, sessInd),sprintf('diagnostic_aamod_tsdiffana%d',i)));
    if (nargin < 4) || ~keepOpen, spm_figure('Close',f(i)); end
end
end
