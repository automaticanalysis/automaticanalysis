function [aap,resp]=aamod_norm_vbm(aap, task, subjInd)
% AAMOD_NORM_VBM Perform's the "old" SPM5/8 segmentation.
%
% [aap,resp] = AAMOD_NORM_VBM(aap, task, subjind)
%
%
% input stream:     structural
%
% output streams:   sn
%                   inv_sn
%                   native_grey
%                   native_white
%                   native_csf
%                   normalised_density_grey
%                   normalised_density_white
%                   normalised_density_csf
%                   normalised_volume_grey
%                   normalised_volume_white
%                   normalised_volume_csf
%

resp='';

switch task
    case 'domain'
        resp='subject';
    case 'description'
        resp='SPM5/8 segmentation for structural images';
    case 'doit'

        
        img = aas_get_files_bystream(aap, subjInd, 'structural');
        
        [pth, nm, ext] = fileparts(img);
               
        estopts.regtype = 'mni';
        estopts.samp = aap.tasklist.currenttask.settings.samp;
        
        out = spm_preproc(img,estopts);
        [sn, isn] = spm_prep2sn(out);        
        matname = fullfile(pth, [nm '_seg_sn.mat']);
        invmatname = fullfile(pth, [nm '_seg_inv_sn.mat']);        
        savefields(matname,sn);
        savefields(invmatname,isn);

        
        % write out GM and WM native + normalised
        writeopts.biascor = 1;
        writeopts.GM  = [1 1 1];	% assume GM(2) means unmod
        writeopts.WM  = [1 1 1];
        writeopts.CSF = [1 1 1];
        
        spm_preproc_write(sn, writeopts);

        % Write normalized images? Only necessary if not using DARTEL etc.
        % for normalization.
        %spm_write_sn(...);


        % Describe outputs
        tiss={'grey','white','csf'};
        for tissind=1:3
            aap = aas_desc_outputs(aap, subjind, sprintf('native_%s', tiss{tissind}), fullfile(pth, sprintf('c%d%s', tissind, [nm ext])));
            aap = aas_desc_outputs(aap, subjind, sprintf('normalised_density_%s', tiss{tissind}), fullfile(pth, sprintf('wc%d%s', tissind, [nm ext])));
            aap = aas_desc_outputs(aap, subjind, sprintf('normalised_volume_%s', tiss{tissind}), fullfile(pth, sprintf('mwc%d%s', tissind, [nm ext])));
        end
        
end
end

%------------------------------------------------------------------------
function savefields(fnam, p)
if length(p)>1
    aas_log(aap, true, 'Can''t save fields.');
end

fn = fieldnames(p);

if numel(fn)==0
    return;
end

for i=1:length(fn),
    feval(@()assignin('caller',fn{i},p.(fn{i})));
end

if str2double(version('-release'))>=14,
    save(fnam,'-V6',fn{:});
else
    save(fnam,fn{:});
end

end
%------------------------------------------------------------------------










