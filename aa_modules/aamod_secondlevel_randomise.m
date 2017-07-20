% AA module - second level statistics
% Only runs if all contrasts present in same order in all subjects at first
% level. If so, makes model with basic t-test for each of contrasts.
% Second-level model from Rik Henson
% Modified for aa by Rhodri Cusack May 2006
%

function [aap,resp]=aamod_secondlevel_randomise(aap,task)

resp='';

switch task
    case 'doit'
        aas_prepare_diagnostic(aap);
        
        nsub = length(aap.acq_details.subjects);
        options = aap.tasklist.currenttask.settings.options;
        
        aas_log(aap,false, sprintf('%d subjects',nsub));
        aas_log(aap, false, sprintf('Will run randomise with options: %s', options));
             
        conFn = cell(1, nsub);
        
        % Get the stream we want to run randomise on
        streams = aap.tasklist.currenttask.inputstreams.stream{:};
        
        % Get data for each subject...
        for subj = 1:nsub            
            % Get the confiles in order...
            conFn{subj} = aas_findstream(aap, streams, subj);
            
            % Get rid of .hdr files if these exist...
            for n = size(conFn{subj}, 1):-1:1
                [pth, nme, ext] = fileparts(deblank(conFn{subj}(n,:)));
                if strcmp(ext, '.hdr')
                    conFn{subj}(n,:) = [];
                end
            end
            
            % We need equal numbers of images across subjects...
            if subj > 1 && size(conFn{subj}, 1) ~= size(conFn{subj - 1}, 1)
                aas_log(aap, true, 'Different number of contrasts across subjects!');
            else
                V = spm_vol(deblank(conFn{subj}(1,:)));
                img_size = V.dim(1) * V.dim(2) * V.dim(3);
            end
        end
        
        aas_log(aap,false,'Merge data for all subjects')
        aapCell = cell(1, size(conFn{1}, 1));
        FSLrandomiseCell = cell(1, size(conFn{1}, 1));
        mergedData = cell(1, size(conFn{1}, 1));
        maskData = cell(1, size(conFn{1}, 1));
        
        pth = aas_getstudypath(aap);
        
        for n = 1:size(conFn{1}, 1)
            aas_log(aap,false,sprintf('Merging data for contrast %d', n))
            
            mergedData{n} = fullfile(pth, sprintf('dataMerged_%04d.nii', n));
            unmergedData = '';
            for subj = 1:nsub
                mask_img([], conFn{subj}(n, :), 0)
                unmergedData = [unmergedData conFn{subj}(n, :) ' '];
            end
            
            FSLmerge = ['fslmerge -t ' mergedData{n} ' ' unmergedData];
            [s, w] = aas_runfslcommand(aap, FSLmerge);
            
            % maskData
            V = spm_vol(mergedData{n});
            Y = spm_read_vols(V);
            V = V(1);
            maskData{n} = fullfile(pth, sprintf('mask_%04d.nii', n));
            V.fname = maskData{n};
            M = all(Y ~= 0, 4);
            spm_write_vol(V, M);
            
            aapCell{n} = aap;
            FSLrandomiseCell{n} = ['randomise -i ' mergedData{n} ' -m ' maskData{n} ' -o ' ... 
                fullfile(pth, sprintf('fslT_%04d', n)) ' -1 ' options];
        end        
        
        switch aap.tasklist.currenttask.settings.parallel
            case {'none', 'serial'}
                for n = 1:length(conFn{1})
                    aas_log(aap,false,sprintf('Running randomise on contrast %d', n))
                    
                    [s, w] = aas_runfslcommand(aap, FSLrandomise{n});
                end
            case 'torque'
                memreq = 64 * img_size * length(conFn);
                timreq = 0.001 * img_size * length(conFn);
                aas_log(aap, false, sprintf('Submitting jobs with %0.2f MB and %0.2f hours', ...
                    memreq/(1024^2), timreq/3600))
                
                [s, w] = qsubcellfun(@aas_runfslcommand, ...
                    aapCell, FSLrandomiseCell, ...
                    'memreq', int32(memreq), ...
                    'timreq', int32(timreq), ...
                    'stack', 1 ...
                    );
        end
        
        aas_log(aap,false,'Clean up merged data')
        for n = 1:size(conFn{1}, 1)
            delete(mergedData{n});
            delete(maskData{n});
        end
        
        %% DECLARE OUTPUTS
        fslt_fns = '';
        fslp_fns = '';
        fslcorrp_fns = '';        
        for n = 1:size(conFn{1}, 1)
            tmp = dir(fullfile(pth, sprintf('fslT_%04d_*_p_tstat*.nii', n)));
            fslp_fns = strvcat(fslp_fns, fullfile(pth, tmp.name));
            
            tmp = dir(fullfile(pth, sprintf('fslT_%04d_*_corrp_tstat*.nii', n)));
            fslcorrp_fns = strvcat(fslcorrp_fns, fullfile(pth, tmp.name));       
            
            tmp = dir(fullfile(pth, sprintf('fslT_%04d_tstat*.nii', n)));
            fslt_fns = strvcat(fslt_fns, fullfile(pth, tmp.name));
            
            %% DIAGNOSTICS (check distribution of T-values in contrasts)
            h = img2hist(fullfile(pth, tmp.name), [], sprintf('con%04d', n));
            saveas(h, fullfile(aap.acq_details.root, 'diagnostics', ...
                [mfilename '_' sprintf('con%04d', n) '.eps']), 'eps');
            try close(h); catch; end
        end
        
        aap=aas_desc_outputs(aap,'secondlevel_fslts', fslt_fns);
        aap=aas_desc_outputs(aap,'secondlevel_fslps', fslp_fns);
        aap=aas_desc_outputs(aap,'secondlevel_fslcorrps', fslcorrp_fns);
        
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end



