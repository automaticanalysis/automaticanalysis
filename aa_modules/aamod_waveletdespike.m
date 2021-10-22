function [aap,resp]=aamod_waveletdespike(aap,task,subj,sess)
resp='';
MAXMEM = 8; % GB

switch task
    case 'report'
        
    case 'doit'
        
        [ ~, SPMtool ] = aas_cache_get(aap,'spm');
        
        [ ~, WDS ] = aas_cache_get(aap,'wds');
        WDS.load;

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
                saveas(f,fullfile(aas_getsesspath(aap,subj,sess),'diagnostic_masking.jpg'));
                
                % save mask
                maskY = maskY > val;
                mask = fullfile(fileparts(mask),'mask.nii');
                nifti_write(mask,maskY,'Mask',maskV);
                thr = 1;
                
            otherwise
                mask = spm_select('FPListRec',SPMtool.toolPath,'brainmask.nii');
        end
        spm_mask(mask,P,thr);
        P = spm_file(P,'prefix','m');
        
        %% Despike
        mem = meminfo;
        args = {'LimitRAM', min([MAXMEM mem.ResFree*0.9])};
        if ~isempty(aas_getsetting(aap,'chainsearch'))
            args = horzcat({'chsearch' aas_getsetting(aap,'chainsearch')},args);
        end
        if ~isempty(aas_getsetting(aap,'threshold'))
            args = horzcat({'threshold' aas_getsetting(aap,'threshold')},args);
        end

        outfilepfx = horzcat(spm_file(P,'path'),repmat(filesep,[size(P,1) 1]),spm_file(P,'basename'));
        WaveletDespike(P,outfilepfx, args{:});
        outputfns = spm_file(outfilepfx,'suffix','_wds.nii');
        if exist(spm_file(outputfns,'ext','nii.gz'),'file'), gunzip(spm_file(outputfns,'ext','nii.gz')); end
        
        %% Describe outputs
        aap=aas_desc_outputs(aap,aap.tasklist.currenttask.domain,[subj sess],streams{1},outputfns);
        
        %% unload
        WDS.unload;
        
    case 'checkrequirements'
        if ~aas_cache_get(aap,'spm'), aas_log(aap,true,'SPM not found'); end
        if ~aas_cache_get(aap,'wds'), aas_log(aap,true,'Wavelet despiking toolbox not found'); end
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end



