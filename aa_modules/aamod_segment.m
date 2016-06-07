% AA module - normalisation using normalise or two pass procedure with segment
% [aap,resp]=aamod_segment(aap, task, subj)
%
% This is an alternative to the traditional aamod_norm_noss. It only
% performs segmentation. If a second-pass of bias correction is wanted, do
% this separately using e.g. aamod_biascorrect_segment8.
function [aap,resp]=aamod_segment(aap, task, subj)
resp='';

switch task
    case 'report'
        
    case 'doit'
                                
        defs = aap.spm.defaults.normalise;
        defs.estimate.weight = '';
        
        
        % Determine the type of image we are expecting (T1, T2, etc.) from
        % the input and output streams specified in the XML file.        
        inStream = aap.tasklist.currenttask.inputstreams.stream{1};
                
        
        structImg = aas_getfiles_bystream(aap, subj, inStream);
        
        % Cheap and cheerful way of ensuring only one file is considered!
        if size(structImg, 1) > 1
            structImg = deblank(structImg(1,:));
            aas_log(aap,0,sprintf('Found more than one structural so using:\n%s',structImg));
        end
        % Get structural directory for this subject
        [structPath, structName, ext] = fileparts(structImg);
        
        
        % Set up normalisation, etc.                    
        estopts.regtype='mni';    % turn on affine again
        V=spm_vol(structImg);
        
        
        % Run normalization
        out = spm_preproc(V,estopts);                            
        [sn,isn]   = spm_prep2sn(out);
                        
        aas_log(aap,false,'Writing out the segmented images')
        % write out GM , WM, CSF native + unmod normalised
        writeopts.biascor = 1;
        writeopts.GM  = [1 1 1];
        writeopts.WM  = [1 1 1];
        writeopts.CSF = [1 1 1];
        writeopts.cleanup = aap.tasklist.currenttask.settings.cleanup;
        spm_preproc_write(sn,writeopts);
            
        SNmat = fullfile(structPath, [structName '_seg_sn.mat']);
        invSNmat = fullfile(structPath, [structName '_seg_inv_sn.mat']);
        savefields(SNmat, sn);
        savefields(invSNmat, isn);
        aap = aas_desc_outputs(aap, subj,'normalisation_seg_sn', SNmat);
        aap = aas_desc_outputs(aap, subj,'normalisation_seg_inv_sn', invSNmat);
        
        % write normalised images        
        spm_write_sn(structImg, SNmat, defs.write);
        
        % describe outputs
        tiss={'grey','white','csf'};
        for tissind=1:3
            aap=aas_desc_outputs(aap, subj, sprintf('normalised_density_%s',tiss{tissind}), fullfile(structPath,sprintf('wc%d%s%s',tissind,structName,ext)));
            aap=aas_desc_outputs(aap, subj, sprintf('normalised_volume_%s',tiss{tissind}), fullfile(structPath,sprintf('mwc%d%s%s',tissind,structName,ext)));
            aap=aas_desc_outputs(aap, subj, sprintf('native_%s',tiss{tissind}), fullfile(structPath,sprintf('c%d%s%s',tissind,structName,ext)));
        end
                        
    case 'checkrequirements'
        
end

%------------------------------------------------------------------------
function savefields(fnam,p)
if length(p)>1, error('Can''t save fields.'); end
fn = fieldnames(p);
if numel(fn)==0, return; end
for subj=1:length(fn),
    feval(@()assignin('caller',fn{subj},p.(fn{subj})));
end
if str2double(version('-release'))>=14,
    save(fnam,'-V6',fn{:});
else
    save(fnam,fn{:});
end

return;
%------------------------------------------------------------------------