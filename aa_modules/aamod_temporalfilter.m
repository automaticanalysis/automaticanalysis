% aa module for temporal filtering using SPM/FieldTrip
% based on Rik Henson's script

function [aap,resp]=aamod_temporalfilter(aap,task,subj,sess)
resp='';

switch task
    case 'report'
        
    case 'doit'
        
        streams=aas_getstreams(aap,'input'); 
        datastream = streams{1};
        headerstream = streams{2};
        
        % obtain input filename(s)
        inputfnames = aas_getfiles_bystream(aap,aap.tasklist.currenttask.domain,[subj sess],datastream);
        pfx = 'f';
        outputfnames = spm_file(inputfnames,'prefix',pfx);
        
        % obtain parameter
        if aas_stream_has_contents(aap,'session',[subj sess],headerstream)
            HDR = load(aas_getfiles_bystream(aap,'session',[subj sess],headerstream));
            TR = HDR.DICOMHEADERS{1}.volumeTR;
        else
            if ~isempty(aap.tasklist.currenttask.settings.TR)
                TR = aap.tasklist.currenttask.settings.TR;
            else
                aas_log(aap,true,'ERROR: TR is not found and not specified!');
            end
        end
        filter  = aap.tasklist.currenttask.settings.filter;
        
        % do stuff (from Rik Henson)
        V = spm_vol(inputfnames);
        Vo = V;
        for i=1:numel(Vo)
            Vo(i).fname = spm_file(V(i).fname,'prefix',pfx);
        end
        Vo = spm_create_vol(Vo);
        
        for i=1:V(1).dim(3)
            dat = spm_read_plane(V,i);
            dat = reshape(dat,prod(V(1).dim(1:2)),[]);
            
            dat = spm_eeg_preproc_filter(aap, filter, dat, 1/TR, i==1);
            
            for j=1:numel(Vo)
                Vo(j) = spm_write_plane(Vo(j),reshape(dat(:,j),V(1).dim(1:2)),i);
            end
        end
        
        % Describe outputs
        aap=aas_desc_outputs(aap,aap.tasklist.currenttask.domain,[subj sess],datastream,outputfnames);
        
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end
end

% ********************************* UTILS *********************************

function Y = spm_read_plane(V,n)

Y = zeros([V(1).dim(1:2) numel(n) numel(V)]);
for i=1:numel(V)
    for j=numel(n)
        Y(:,:,j,i) = spm_slice_vol(V(i),spm_matrix([0 0 n(j)]),V(i).dim(1:2),0);
    end
end
if numel(n)==1, Y = Y(:,:,1,:); end
end

function dat = spm_eeg_preproc_filter(aap, filter, dat, Fs, verbose)

Fp   = [filter.LowCutoffFreqency filter.HighCutoffFreqency];
type = filter.type;
N    = filter.order;
dir  = filter.direction;

switch numel(Fp)
    case 2
        if filter.bandstop
            filter_func = @ft_preproc_bandstopfilter;
            if verbose, aas_log(aap,false,'INFO: Two cutoff frequencies provided --> bandstop will be applied'); end
        else
            filter_func = @ft_preproc_bandpassfilter;
            if verbose, aas_log(aap,false,'INFO: Two cutoff frequencies provided --> bandpass will be applied'); end
        end
    case 1
        if isempty(filter.HighCutoffFreqency) % highpass
            filter_func = @ft_preproc_highpassfilter;
            if verbose, aas_log(aap,false,'INFO: One cutoff frequency (LowCutoffFrequency) provided --> highpass will be applied'); end
        else % lowpass
            filter_func = @ft_preproc_lowpassfilter;
            if verbose, aas_log(aap,false,'INFO: One cutoff frequency (HighCutoffFreqency) provided --> lowpass will be applied'); end
        end
    otherwise
        aas_log(aap,true,'ERROR: No cutoff frequency provided!');
end

dat = filter_func(dat,Fs,Fp,N,type,dir);

end