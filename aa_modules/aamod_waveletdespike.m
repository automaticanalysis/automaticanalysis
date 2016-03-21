function [aap,resp]=aamod_waveletdespike(aap,task,subj,sess)
resp='';
MAXMEM = 8; % GB

switch task
    case 'report'
        
    case 'doit'
        streams=aap.tasklist.currenttask.inputstreams.stream;
        
        % Images to smooth (1st stream)
        P = aas_getfiles_bystream(aap,aap.tasklist.currenttask.domain,[subj sess],streams{1});
        
        %% masking
        thr = aap.tasklist.currenttask.settings.maskingthreshold;
        switch thr
            case 1
                mask = aas_getfiles_bystream_multilevel(aap,aap.tasklist.currenttask.domain,[subj,sess],'brainmask');
                
            case 0
                % read image
                mask = [P(1,:),',1'];
                maskV = spm_vol(mask);
                maskY = spm_read_vols(maskV);
                
                % estimate threshold
                [n, x] = hist(maskY(:),100);
                dn = abs(diff(n));
                ind = find(dn<10,1,'first');
                val = x(ind+1);
                
                % diag
                f = figure; hold on;
                plot(x,n); plot([val val],[0 max(n)],'r');
                title('Intensity distribution with masking threshold');
                saveas(f,fullfile(aas_getsesspath(aap,subj,sess),'diagnostics_masking.jpg'));
                
                % save mask
                maskY = maskY > val;
                mask = fullfile(fileparts(mask),'mask.nii');
                nifti_write(mask,maskY,'Mask',maskV);
                thr = 1;
                
            otherwise
                mask = spm_select('FPListRec',aap.directory_conventions.spmdir,'brainmask.nii');
        end
        spm_mask(mask,P,thr);
        P = spm_file(P,'prefix','m');
        
        %% Despike
        if isfield(aap.tasklist.currenttask.settings,'chainsearch')
            chsearch = aap.tasklist.currenttask.settings.chainsearch;
        else
            chsearch = 'moderate';
        end
        mem = meminfo;
        outfilepfx = horzcat(spm_file(P,'path'),repmat(filesep,[size(P,1) 1]),spm_file(P,'basename'));
        WaveletDespike(P,outfilepfx, ...
            'chsearch', chsearch,...
            'LimitRAM', min([MAXMEM mem.MemFree*0.9]));
        outputfns = spm_file(outfilepfx,'suffix','_wds.nii');
        
        %% Describe outputs
        aap=aas_desc_outputs(aap,aap.tasklist.currenttask.domain,[subj sess],streams{1},outputfns);
        
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end;



