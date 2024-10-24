function [ aap,resp ] = aamod_secondlevel_randomise(aap,task)

% AA module - nonparametric secondlevel statistics using FSL randomise
%
% Only runs if all contrasts present in same order in all subjects at first
% level. If so, makes model with basic t-test for each of contrasts.
%
% NEW FILE NAMING CONVENTION
%
%  FSL randomise doesn't exactly have the most transparent output filenaming  
%  Ergo, we attempt to convert them to something more readable  
%
%  In brief: fslT_000C_tstat1.nii => <cname>_<identifier>_uncorrected.nii  
%  and the thresholded version = <cname>_<identifier>_<pmap>_corrected.nii  
%  where C is a contrast number, cname is the corresponding contrast name,  
%  "identifier" is an optional identifier (e.g., "young" or "controls")  
%  and pmap is the threhsolding option  
%
%  Additionally, if the glm_output option was passed to aamod_secondlevel_randomise  
%  the generated cope file is reanamed to <cname>_<identifier>_cope.nii. The sigmasqr,  
%  pe, and var files are deleted (UPDATE  we prolly can keep, since the new directory
% structure is less clusttered)  
% 
%  if multiple tstat files are present (e.g., *_tstat1.nii, *_tstat2.nii  
%  as generated in a 2-group ttest) the new names include the tstat number  
%
% CHANGE HISTORY
%
% 10/2024 [MSJ] swap "sleep" for "pause" and { } for "extractfield"
% 08/2023 [MSJ] skip F contrasts
% 06/2021 [MSJ] attempt to implement sensible FSL file naming
% 06/2019 [MSJ] Modified to run one or two sample ttest w/ covariates
%
% Modified for aa by Rhodri Cusack May 2006
% Second-level model from Rik Henson

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
        master_subject_list = { aap.acq_details.subjects.subjname };
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
        
        FSLCommandString = cell(1, size(conFn{1}, 1));
        mergedData = cell(1, size(conFn{1}, 1));
        maskData = cell(1, size(conFn{1}, 1));
        
        pth = aas_getstudypath(aap);
        
        % sort files by contrast, each in a separate folder
        % named according to the contrast to better organize results

        % need an SPM to get the contrast names, which we assume are
        % the same for all subjects, so we can just use first subject

        SPM = load(aas_getfiles_bystream(aap, 1, 'firstlevel_spm'));
        SPM = SPM.SPM;
        dirnames = { SPM.xCon(:).name };       

        % aa weirdness: the stream firstlevel_cons only contains T
        % contrasts (because it only saves con_xxxx.nii and
        % SPM writes F contrasts to ess_xxxx.nii). Currently,
        % we simply ignore any F contrasts

        fcheck = strcmp({SPM.xCon.STAT},'F');

        if any(fcheck)
            aas_log(aap, false, sprintf('Ignoring F contrasts'));
            dirnames(fcheck) = [];
        end

        % take out characters that don't go well in filenames...
        dirnames = cellfun(@(x) char(regexp(x,'[a-zA-Z0-9_-]','match'))', dirnames, 'UniformOutput',false);

        % output file extension may be .nii or nii.gz
        % BUGWATCH: need to deal with nii.gz in aamod_secondlevel_randomise_threshold

        if (strcmp(aap.directory_conventions.fsloutputtype,'NIFTI') ~= 1)
            fsl_file_extension = 'nii.gz';
        else
            fsl_file_extension = 'nii';
        end        
        
        for cindex = 1:numel(dirnames)
               
            % recall we set pwd to module directory on entry...
            dirname = dirnames{cindex};
            aas_makedir(aap,dirname); % harmless if dirname exists
            pause(1);   % make sure filesystem catches up
            
            aas_log(aap,false,sprintf('Merging data for contrast %d', cindex))
            
            mergedData{cindex} = fullfile(pth, sprintf('%s/mergedata.nii', dirname));
          
            unmergedData = '';
            for sindex = 1:nsub
                mask_img([], conFn{sindex}(cindex, :), 0)
                unmergedData = [unmergedData conFn{sindex}(cindex, :) ' '];
            end
            
            FSLmerge = ['fslmerge -t ' mergedData{cindex} ' ' unmergedData];
            aas_runfslcommand(aap, FSLmerge);
            
            % make a mask on the fly
            
            V = spm_vol(mergedData{cindex});
            Y = spm_read_vols(V);
            V = V(1);
            maskData{cindex} = fullfile(pth, sprintf('%s/mask.nii', dirname));
            V.fname = maskData{cindex};
            M = all(Y ~= 0, 4);
            spm_write_vol(V, M);
            
            % generate a Randomise command for the test we're doing
            
            % if we're doing backgrounding, terminate the randomize command with a bg token

            if (aap.tasklist.currenttask.settings.runFSLinbackground)
                options = [ aap.tasklist.currenttask.settings.options ' &'];
            else
                options = aap.tasklist.currenttask.settings.options;
            end

            % regardless of bgoption, DON'T run the final contrast in the background
            % (we need one job to block (hopefully all other jobs will complete)

            if (cindex == numel(dirnames))
                options = aap.tasklist.currenttask.settings.options;
            end

            % we need a unique separator for a_2_r_threshold to find
            % use double underscore + identifier + double underscore
            % if user didn't specify an identifier, use XXX
            % note we only specify the first trailing underscore here
            % -- FSL will add the second

            if (isempty(aap.tasklist.currenttask.settings.identifier))
                outfile_prefix = [ dirname '__XXX_'];
            else
                outfile_prefix = [ dirname '__' aap.tasklist.currenttask.settings.identifier '_'];
            end
               
            if ((number_of_testgroups == 1) && isempty(covariates))
                
                % special simple case: 1 sample ttest
                % reference: randomise -i dataMerged_0001.nii -o ttest1 -1 -T -v 5 -n 5000
                
                FSLCommandString{cindex} = ['randomise -i ' mergedData{cindex} ' -m ' maskData{cindex} ' -o ' ... 
                    fullfile(pth, sprintf('%s/%s',dirname,outfile_prefix)) ' -1 ' options];

            else

                % reference: randomise -i dataMerged_0001.nii -o ttest2 -d design.mat -t design.con -c 3.1 -x -n 1000
                
                SID_1 = aap.tasklist.currenttask.settings.group_one_subjectIDs;
                SID_2 = aap.tasklist.currenttask.settings.group_two_subjectIDs;
  
                % this will create design.mat and design.con
                % could move this outside con loop (design.mat and
                % design.con are the same for all contrasts...) 
                 
                create_auxillary_fsl_files(aap, master_subject_list, SID_1, SID_2, covariates);
                
                FSLCommandString{cindex} = ['randomise -i ' mergedData{cindex} ' -m ' maskData{cindex} ' -o ' ... 
                    fullfile(pth, sprintf('%s/%s',dirname,outfile_prefix)) ' -d ' fullfile(pth,'design.mat') ' -t ' fullfile(pth,'design.con') ' ' options];

            end
            
            
        end   % loop over covariates        
        
        % now run FSL
                
        for cindex = 1:size(conFn{1}, 1)
            aas_log(aap,false,sprintf('Running randomise on contrast %d', cindex))
            aas_runfslcommand(aap, FSLCommandString{cindex});
            pause(1);  % feels necessary to let bg process launch...
        end
              
        % block until all the contrast directories are populated

        for cindex = 1:numel(dirnames)

            dirname = dirnames{cindex};

            temp = dir(fullfile(pth, sprintf('%s/*tstat*.%s', dirname, fsl_file_extension)));

            while (numel(temp) == 0)
                pause(60);
                temp = dir(fullfile(pth, sprintf('%s/*tstat*.%s', dirname, fsl_file_extension)));
             end

        end
        
        % one final block here in case filesystem needs to catch-up
        pause(60);

        % clean up temp files; Note we don't delete design.mat and 
        % design.con in two-group test -- leave for verification
        
        for cindex = 1:size(conFn{1}, 1)
            delete(mergedData{cindex});
            delete(maskData{cindex});
        end
        
        % desc the outputs
        
        fslt_fns = '';
        fslcorrp_fns = '';
        fslcope_fns = '';
        
        pindex = 1;
        tindex = 1;
        gindex = 1;

        for cindex = 1:numel(dirnames)
                        
            dirname = dirnames{cindex};
                                  
            if (isempty(aap.tasklist.currenttask.settings.identifier))
                outfile_prefix = [ dirname '__XXX_'];
            else
                outfile_prefix = [ dirname '__' aap.tasklist.currenttask.settings.identifier '_'];
            end
          
            % uncorrected tmap (./conname/prefix_tstat.extension)
            
            temp = dir(fullfile(pth, sprintf('%s/%s_tstat*.%s', dirname, outfile_prefix, fsl_file_extension)));

            if (isempty(fslt_fns)); fslt_fns = cell(size(conFn{1},1)*size(temp,1),1); end
            
            for findex = 1:size(temp,1)
                fslt_fns{tindex} = [temp(findex).folder '/' temp(findex).name];
                tindex = tindex + 1;
            end
            
            % various p-correction files based on options (./conname/prefix_*corrp*.nii[.gz])
            
            temp = dir(fullfile(pth, sprintf('%s/%s*corrp*.%s', dirname, outfile_prefix, fsl_file_extension)));
            
            if (isempty(fslcorrp_fns)); fslcorrp_fns = cell(size(conFn{1},1)*size(temp,1),1); end
            
            for findex = 1:size(temp,1)
                fslcorrp_fns{pindex} = [temp(findex).folder '/' temp(findex).name];
                pindex = pindex + 1;
            end
            
           % cope files if --glm_output was used (./conname/prefix_*glm*.nii[.gz])
            
            temp = dir(fullfile(pth, sprintf('%s/%s*glm*.%s', dirname, outfile_prefix, fsl_file_extension)));
            
            if (isempty(fslcope_fns)); fslcope_fns = cell(size(conFn{1},1)*size(temp,1),1); end
            
            for findex = 1:size(temp,1)
                fslcope_fns{gindex} = [temp(findex).folder '/' temp(findex).name];
                gindex = gindex + 1;
            end
            
        end
        
        aap=aas_desc_outputs(aap,'secondlevel_fslts', fslt_fns);
        aap=aas_desc_outputs(aap,'secondlevel_fslcorrps', fslcorrp_fns);
        aap=aas_desc_outputs(aap,'secondlevel_fslcopes', fslcope_fns);
     
        % document subjects and covariates for QA

        aas_prepare_diagnostic(aap);
         
        writetable(cell2table(master_subject_list'), fullfile(aas_getstudypath(aap),'diagnostics','diagnostic_SUBJECTS.txt'), 'WriteVariableNames', false);

        if (~isempty(covariates))
            writetable(array2table(covariates), fullfile(aas_getstudypath(aap), 'diagnostics', 'diagnostic_COVARIATES.txt')); 
        end   
                
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

    % the only time this if branch should be entered is one group + one
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
