% AA module - realignment and unwarp
% [aap,resp]=aamod_realignunwarp(aap,task,i)
% Realignment using SPM5
% i=subject num
% Rhodri Cusack MRC CBU 2004-6 
% @@@ THIS IS NOT YET TRANSFORMED TO AA4 @@@

function [aap,resp]=aamod_realignunwarp(aap,task,i)

resp='';

switch task
    case 'domain'
        resp='subject';
    case 'description'
        resp='SPM5 realign and unwarp';
    case 'summary'
        resp='Done SPM5 realign and unwarp\n';
    case 'report'
        mvmean=[];
        mvmax=[];
        mvstd=[];
        mvall=[];
        nsess=length(aap.acq_details.sessions);

        qq=[];
        for j=1:nsess
            im1fn=aas_getimages(aap,i,j,aap.tasklist.currenttask.epiprefix,aap.acq_details.numdummies,1+aap.acq_details.numdummies);
            im1V=spm_vol(im1fn);
            qq(j,:)     = spm_imatrix(im1V.mat);
            rpfn=spm_select('List',aas_getsesspath(aap,i,j),'^rp.*txt');
            mv=spm_load(fullfile(aas_getsesspath(aap,i,j),rpfn));
            mv=mv+repmat(qq(j,1:6)-qq(1,1:6),[size(mv,1) 1]);
            mv(:,4:6)=mv(:,4:6)*180/pi; % convert to degrees!
            mvmean(j,:)=mean(mv);
            mvmax(j,:)=max(mv);
            mvstd(j,:)=std(mv);
            mvall=[mvall;mv];
        end;

        
        aap.report.html=strcat(aap.report.html,'<h3>Movement maximums</h3>');
        aap.report.html=strcat(aap.report.html,'<table cellspacing="10">');
        aap.report.html=strcat(aap.report.html,sprintf('<tr><td align="right">Sess</td><td align="right">x</td><td align="right">y</td><td align="right">z</td><td align="right">rotx</td><td align="right">roty</td><td align="right">rotz</td></tr>',j));
        for j=1:nsess
        aap.report.html=strcat(aap.report.html,sprintf('<tr><td align="right">%d</td>',j));
        aap.report.html=strcat(aap.report.html,sprintf('<td align="right">%8.3f</td>',mvmax(j,:)));
        aap.report.html=strcat(aap.report.html,sprintf('</tr>',j));
        end;
        aap.report.html=strcat(aap.report.html,'</table>');
        
        varcomp=mean((std(mvall).^2)./(mean(mvstd.^2)));
        aap.report.html=strcat(aap.report.html,'<h3>All variance vs. within session variance</h3><table><tr>');
        aap.report.html=strcat(aap.report.html,sprintf('<td>%8.3f</td>',varcomp));
        aap.report.html=strcat(aap.report.html,'</tr></table>');
        
        
        
        
        aap=aas_report_addimage(aap,fullfile(aas_getsubjpath(aap,i),'diagnostic_aamod_realign.jpg'));     

    case 'doit'
        
        % Load up template job description
        load spmjob_realignunwarp.mat
        
        jobs{1}.spatial{1}.realignunwarp.data=[];
        for j = aap.acq_details.selected_sessions % 
            clear mydata;
            % get files in this directory
            P = aas_getimages(aap,i,j,aap.tasklist.currenttask.epiprefix);
            mydata{j} = P;
            jobs{1}.spatial{1}.realignunwarp.data(j).scans=mydata;
            jobs{1}.spatial{1}.realignunwarp.data(j).pmscan={};
        end

        
        spm_jobman('run',jobs);
        
        % Save graphical output
        try figure(spm_figure('FindWin', 'Graphics')); catch; figure(1); end;
        print('-djpeg','-r75',fullfile(aas_getsubjpath(aap,i),'diagnostic_aamod_realignunwarp'));
                
case 'checkrequirements'
    
otherwise
    aas_log(aap,1,sprintf('Unknown task %s',task));
end;
