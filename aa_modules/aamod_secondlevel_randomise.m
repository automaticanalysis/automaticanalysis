% AA module - nonparametric secondlevel statistics using FSL randomise
%
% Only runs if all contrasts present in same order in all subjects at first
% level. If so, makes model with basic t-test for each of contrasts.
%
% Modified for aa by Rhodri Cusack May 2006
% Second-level model from Rik Henson
%
% CHANGE HISTORY
%
% 06/2019 [MSJ] Modified to run one or two sample ttest w/ covariates
%

function [ aap,resp ] = aamod_secondlevel_randomise(aap,task)

resp='';

switch task
    
    case 'doit'
        
        savedir = pwd;
        cd(aas_getstudypath(aap));
        
        % figure out what kind of test we're doing...
        
        number_of_testgroups = 1;
        
        if (~isempty(aap.tasklist.currenttask.settings.group_two_subjectIDs))
            number_of_testgroups = 2;
        end
            
        % sanity checks
        
        if (number_of_testgroups == 2 && isempty(aap.tasklist.currenttask.settings.group_one_subjectIDs))
             aas_log(aap, true, sprintf('You must explicitly specify both subject groups in a two-sample ttest. Exiting...'));
        end
        
        if (~isempty(intersect(aap.tasklist.currenttask.settings.group_one_subjectIDs,aap.tasklist.currenttask.settings.group_two_subjectIDs)))
              aas_log(aap, true, sprintf('At least one subject ID is on both group lists. Exiting...'));
        end  
        
        % prepare the list of subjects to be included
        
        master_subject_list = {};
                
        if ~isempty(aap.tasklist.currenttask.settings.group_one_subjectIDs)
            temp = aap.tasklist.currenttask.settings.group_one_subjectIDs;
            [ r,c ] = size(temp);
            if (r>c); temp = temp'; end
            master_subject_list = temp;
            if ~isempty(aap.tasklist.currenttask.settings.group_two_subjectIDs)
                temp = aap.tasklist.currenttask.settings.group_two_subjectIDs;
                [ r,c ] = size(temp);
                if (r>c); temp = temp'; end
                master_subject_list = { master_subject_list{:}, temp{:} };
            end   
        else
            master_subject_list = extractfield(aap.acq_details.subjects,'subjname');
        end
        
        nsub = numel(master_subject_list);
        
        if (~isempty(aap.tasklist.currenttask.settings.covariates))
            covariates = process_covariate_option(aap, aap.tasklist.currenttask.settings.covariates, aap.tasklist.currenttask.settings.demean_covariates);
            if (size(covariates,1) ~= nsub)
                % saving throw -- try transposing covariate array (user may have entered row-wise instead of column-wise)
                covariates = transpose(covariates);
                if (size(covariates,1) ~= nsub)
                    aas_log(aap, true, sprintf('There must be one covariate entry for each subject. Exiting...'));
                end
            end
        else
            covariates = [];
        end

        options = aap.tasklist.currenttask.settings.options;
        
        aas_log(aap,false, sprintf('Running FSL Randomise on %d subjects with options: %s', nsub, options));
                      
        % Get con files for each subject
        
        conFn = cell(1, nsub);
        
        % careful -- aas_getfiles_bystream takes a subject INDEX (in the
        % order the subjects were added) so we have to loop over all of
        % them and check the name subfield against the subject names
        
        for sindex = 1:length(aap.acq_details.subjects)
                    
            master_index = find(strcmp(master_subject_list,aap.acq_details.subjects(sindex).subjname));
            
            if (~isempty(master_index) && master_index>0)
                 conFn{master_index} = aas_getfiles_bystream(aap, sindex, 'firstlevel_cons');
                 if (isempty(conFn{master_index}))
                    badname = aap.acq_details.subjects(sindex).subjname;
                    aas_log(aap, true, sprintf('Empty contrast returned for subject %s (sindex %d). Check group_one_subjectIDs and/or group_two_subjectIDs for invalid entries.', badname, sindex));
                 end
            end
           
        end

        % sanity check -- look for holes in conFn 
        % this can happen if there's bad subject IDs in group_xxx_subjectIDs
        % and checking now is better than mysteriously crashing later
        
        for sindex = 1:numel(conFn)
           if (isempty(conFn{sindex}))
                aas_log(aap, true, sprintf('Empty contrast found (sindex %d). Check group_one_subjectIDs and/or group_two_subjectIDs for invalid entries', sindex));
           end
        end
          
              
        % FSL requires we merge all the data into one big 4D file
        % -- the same file can be used for a one sample or two
        % sample ttest etc; we simply pass a different design matrix
        
        FSLrandomiseCell = cell(1, size(conFn{1}, 1));
        mergedData = cell(1, size(conFn{1}, 1));
        maskData = cell(1, size(conFn{1}, 1));
        
        pth = aas_getstudypath(aap);
            
        for n = 1:size(conFn{1}, 1)
            
            aas_log(aap,false,sprintf('Merging data for contrast %d', n))
            
            mergedData{n} = fullfile(pth, sprintf('dataMerged_%04d.nii', n));
            unmergedData = '';
            for sindex = 1:nsub
                mask_img([], conFn{sindex}(n, :), 0)
                unmergedData = [unmergedData conFn{sindex}(n, :) ' '];
            end
            
            FSLmerge = ['fslmerge -t ' mergedData{n} ' ' unmergedData];
            aas_runfslcommand(aap, FSLmerge);
            
            % make a mask on the fly
            
            V = spm_vol(mergedData{n});
            Y = spm_read_vols(V);
            V = V(1);
            maskData{n} = fullfile(pth, sprintf('mask_%04d.nii', n));
            V.fname = maskData{n};
            M = all(Y ~= 0, 4);
            spm_write_vol(V, M);
            
            % generate a Randomise command for the test we're doing
                
            if ((number_of_testgroups == 1) && isempty(covariates))
                
                % special simple case: 1 sample ttest
                % reference: randomise -i dataMerged_0001.nii -o ttest1 -1 -T -v 5
                
                FSLrandomiseCell{n} = ['randomise -i ' mergedData{n} ' -m ' maskData{n} ' -o ' ... 
                    fullfile(pth, sprintf('fslT_%04d', n)) ' -1 ' options];

            else

                % reference: randomise -i dataMerged_0001.nii -o ttest2 -d design.mat -t design.con -c 3.1 -x -n 1000
                
                SID_1 = aap.tasklist.currenttask.settings.group_one_subjectIDs;
                SID_2 = aap.tasklist.currenttask.settings.group_two_subjectIDs;
  
                % this will create design.mat and design.con
                % could move this outside con loop (design.mat and
                % design.con are the same for all contrasts...)
                
                create_auxillary_fsl_files(aap, master_subject_list, SID_1, SID_2, covariates);
                 
                FSLrandomiseCell{n} = ['randomise -i ' mergedData{n} ' -m ' maskData{n} ' -o ' ... 
                    fullfile(pth, sprintf('fslT_%04d', n)) ' -d ' fullfile(pth,'design.mat') ' -t ' fullfile(pth,'design.con') ' ' options];

            end
            
        end     % loop over covariates        
        
        % now run FSL
                
        for n = 1:size(conFn{1}, 1)
            aas_log(aap,false,sprintf('Running randomise on contrast %d', n))
            aas_runfslcommand(aap, FSLrandomiseCell{n});
        end

        % clean up temp files; Note we don't delete design.mat and 
        % design.con in two-group test -- leave for verification
        
        for n = 1:size(conFn{1}, 1)
            delete(mergedData{n});
            delete(maskData{n});
        end
        
        % desc the outputs and save some diagnostic images
        
        fslt_fns = '';
        fslcorrp_fns = '';
        
        pindex = 1;
        tindex = 1;
        
        aas_prepare_diagnostic(aap);

        for contrast = 1:size(conFn{1}, 1)
                        
            tmp = dir(fullfile(pth, sprintf('fslT_%04d_*_corrp_tstat*.nii', contrast)));
            
            if (isempty(fslcorrp_fns)); fslcorrp_fns = cell(size(conFn{1},1)*size(tmp,1),1); end
            
            for dindex = 1:size(tmp,1)
                fslcorrp_fns{pindex} = [tmp(dindex).folder '/' tmp(dindex).name];
                pindex = pindex + 1;
            end
            
            tmp = dir(fullfile(pth, sprintf('fslT_%04d_tstat*.nii', contrast)));

            if (isempty(fslt_fns)); fslt_fns = cell(size(conFn{1},1)*size(tmp,1),1); end
            
            for dindex = 1:size(tmp,1)
                fslt_fns{tindex} = [tmp(dindex).folder '/' tmp(dindex).name];
                tindex = tindex + 1;
            end
            
        end
        
        % document subjects and covariates for QA
        
        writetable(cell2table(master_subject_list'), fullfile(aas_getstudypath(aap),'diagnostic_SUBJECTS.txt'), 'WriteVariableNames', false);

        if (~isempty(covariates))
            writetable(array2table(covariates), fullfile(aas_getstudypath(aap),'diagnostic_COVARIATES.txt')); 
        end

        aap=aas_desc_outputs(aap,'secondlevel_fslts', fslt_fns);
        aap=aas_desc_outputs(aap,'secondlevel_fslcorrps', fslcorrp_fns);
        
        % restore pwd and we're done!

        cd(savedir);
        
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
	end
        
end




%--------------------------------------------------------------------------------------------------------------------
function create_auxillary_fsl_files(aap, master_subject_list, group_one_subjectIDs, group_two_subjectIDs, covariates)
%--------------------------------------------------------------------------------------------------------------------

design_matrix = [];
temp = [];

if (isempty(group_two_subjectIDs))
    
    % one sample ttest with optional covariate(s)
    
    for index = 1:numel(master_subject_list)
        
        if (~isempty(covariates)); temp = covariates(index,:); end
        
        design_matrix = [ design_matrix ; [ 1 temp ] ];
                
    end

    % this only time this if branch should be entered is one group + one
    % covariate (one group no covariates uses simplified fsl call (see line
    % ~139 and two groups uses branch below). As such, an eye(2) contrast
    % isn't useful here -- if the covariate is nuisance, we just want [1 0]
    % [1 0 ...; -1 0 ...] and if the covarate is of interest (i.e. whole brain
    % corr) we want [0 1; 0 -1]

    contrast_matrix = [ 0 1; 0 -1 ];

else
    
    % two sample ttest with optional covariate(s)

    for index = 1:numel(master_subject_list)
        
       if (~isempty(covariates)); temp = covariates(index,:); end
       
         if (any(strcmp(group_one_subjectIDs, master_subject_list{index})))
            design_matrix = [ design_matrix ; [ 1 0 temp ] ];
        end
        
        if (any(strcmp(group_two_subjectIDs, master_subject_list{index})))
            design_matrix = [ design_matrix ; [ 0 1 temp ] ];
        end
        
    end
        
    contrast_matrix = eye(2+size(covariates,2)); % empty "covariates" ok here
    contrast_matrix(1,2) = -1;
    contrast_matrix(2,1) = -1;

end

    
dlmwrite('design.txt', design_matrix, 'delimiter', ' ');
dlmwrite('contrast.txt', contrast_matrix, 'delimiter', ' ');

% runfslcommand will halt on error -- no need to check return flag

aas_runfslcommand(aap, 'Text2Vest design.txt design.mat');
aas_runfslcommand(aap, 'Text2Vest contrast.txt design.con');

end



%-------------------------------------------------------------------------------------------
function covariates = process_covariate_option(aap, specifier, do_mean_centering)
%-------------------------------------------------------------------------------------------

covariates = [];

if (isempty(specifier)); return; end

if (isnumeric(specifier))
    covariates = specifier;
    if (do_mean_centering)
        covariates = covariates - mean(covariates); 
    end
    return;
end

if (ischar(specifier))
    if (~exist(specifier,'file'))
        aas_log(aap, true, sprintf('Could not read covariate file %s. Exiting...', specifier)); 
    end
    % 'importdata' expects txt extension...
    [ ~,~,e ] = fileparts(specifier);
    if (~strcmp(e,'.txt'))
        aas_log(aap, true, 'Covariate file must have ''.txt'' extension. Exiting...'); 
    end
    temp = importdata(specifier);
    if (isstruct(temp))
        covariates = temp.data;
    else
        covariates = temp;
    end
    if (do_mean_centering)
        covariates = covariates - mean(covariates);
    end
    return;
end

aas_log(aap, true, sprintf('Unrecoginzed covariate specifier. Exiting... %s', specifier));

end
