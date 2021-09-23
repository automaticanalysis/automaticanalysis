% aa module - Mask images from a segmented structural
% 1) Take segmentation images
% 2) Threshold them at one of the 3 different levels (zero level, strict [e.g. 99%],
% or an exclusive [highest tissue probability wins])

function [aap,resp] = aamod_mask_fromsegment(aap,task,subj)

resp = '';

switch task
    case 'doit'
        
        inStreams = aas_getstreams(aap,'input');
        refimg = ''; ind = 1:numel(inStreams);
        if aas_stream_has_contents(aap,inStreams{1})
            refimg = aas_getfiles_bystream(aap,subj,inStreams{1});
        end
        ind(1) = [];
        segimg = arrayfun(@(x) aas_getfiles_bystream(aap,subj,inStreams{x}),ind,'UniformOutput',false);        
        
        %% 1) RESLICE
        if ~isempty(refimg)
            % Get defaults
            resFlags = aap.spm.defaults.coreg.write;
            resFlags.which = 1;
            resFlags.mean = 0;
            resFlags.mask = 1;
            
            spm_reslice([{refimg},segimg],resFlags);
            segimg = spm_file(segimg,'prefix',resFlags.prefix);
        end
        
        %% 2) THRESHOLD
        
        %% IN PROGRESS
        outstream = '';
        
        % Load the correct resliced file!
        V = spm_vol(segimg); V = cell2mat(V);
        Y = spm_read_vols(V);
        
        switch aap.tasklist.currenttask.settings.threshold
            case 'zero'
                %% A) Zero thresholding is easy...
                for a = 1:numel(segimg)
                    Y(:,:,:,a) = Y(:,:,:,a) > 0;
                    
                    V(a).fname = spm_file(segimg{a},'prefix','Z_');
                    outstream = strvcat(outstream, V(a).fname); % Save to stream...
                    spm_write_vol(V(a), Y(:,:,:,a));
                    
                    aas_log(aap,false,sprintf('Zero thresholded image %s sums up to %d vox', ...
                        V(a).fname, sum(Y(:,:,:,a),'all')))
                end
            case 'exclusive'
                %% C) Exclusive thresholding of masks:
                % Any particular voxel has greatest chance of being...
                maxeY = max(Y,[],4);
                
                for a = 1:numel(segimg)
                    % We want to check where the segmentation has the
                    % greatest values
                    Y(:,:,:,a) = (Y(:,:,:,a) == maxeY) & (Y(:,:,:,a) > 0.01);
                    
                    V(a).fname = spm_file(segimg{a},'prefix','E_');
                    outstream = strvcat(outstream, V(a).fname); % Save to stream...
                    spm_write_vol(V(a), Y(:,:,:,a));
                    
                    aas_log(aap,false,sprintf('Exclusive thresholded image %s sums up to %d vox', ...
                        V(a).fname, sum(Y(:,:,:,a),'all')))
                end
            otherwise
                %% B) Specific thresholding of each mask
                thr = aas_getsetting(aap,'threshold');
                if numel(thr) == 1, thr(1:3) = thr; end
                for a = 1:numel(segimg)
                    maxY = max(Y(:,:,:,a),[],'all');
                    Y(:,:,:,a) = Y(:,:,:,a) > thr(a);
                    
                    if numel(Y(:,:,:,a) > 0) == 0 % check for bad thesholding
                        aas_log(aap, true, sprintf('ERROR: No voxels above the threshold mask (%f) [max: %f]', thr(a), maxY))
                    end
                    
                    V(a).fname = spm_file(segimg{a},'prefix','S_');
                    outstream = strvcat(outstream, V(a).fname); % Save to stream...
                    spm_write_vol(V(a), Y(:,:,:,a));
                    
                    aas_log(aap,false,sprintf('Strict thresholded image %s sums up to %d vox', ...
                        V(a).fname, sum(Y(:,:,:,a),'all')))
                end
        end
        
        %% OUTPUT
        outStreams = aas_getstreams(aap,'output');
        for a = 1:numel(inStreams)-1 % not for reference
            aap = aas_desc_outputs(aap,subj,outStreams{a},outstream(a,:));
        end
        aap = aas_desc_outputs(aap,subj,'segmasks',outstream);
    case 'checkrequirements'
        in = aas_getstreams(aap,'input'); in(1) = []; % not for reference
        [stagename, index] = strtok_ptrn(aap.tasklist.currenttask.name,'_0');
        stageindex = sscanf(index,'_%05d');
        out = aap.tasksettings.(stagename)(stageindex).outputstreams.stream; if ~iscell(out), out = {out}; end
        for s = 1:numel(in)
            instream = textscan(in{s},'%s','delimiter','.'); instream = instream{1}{end};
            if ~strcmp(out{s},[instream '_mask'])                
                aap = aas_renamestream(aap,aap.tasklist.currenttask.name,out{s},[instream '_mask'],'output');
                aas_log(aap,false,['INFO: ' aap.tasklist.currenttask.name ' output stream: ''' [instream '_mask'] '''']);
            end
        end
end
