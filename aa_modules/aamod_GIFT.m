% Run a group ICA analysis with 'gift'

function [aap, resp] = aamod_groupICA(aap, task, subj)

resp='';

switch task
    case 'domain'
        resp='study';   % this module needs to be run once per study
        
    case 'description'
        resp='Group ICA analysis';
        
    case 'summary'
        subjpath = aas_getsubjpath(subj);
        resp = sprintf('Gift model %s\n', subjpath);
        
    case 'report'
        
    case 'doit'

        settings = aap.tasklist.currenttask.settings;
        
        defaultParmFile = settings.parameterFile;
        if ~exist(defaultParmFile), error('Can''t find %s!', defaultParmFile); end
        load(defaultParmFile);
           
        subjPath = aas_getsubjpath(aap, subj);
        parmFile = '_ica_parameter_info.mat';
        
        sesInfo.userInput.pwd = subjPath;
        sesInfo.userInput.param_file = fullfile(subjPath, parmFile);
        
        nscans = [];
        files = struct('name', {});       
        for ses = 1 : length(aap.acq_details.selected_sessions)
            files(end+1) = struct('name', aas_getimages_bystream(aap, subj, ses, 'epi'));
            nscans(end+1) = size(files(ses).name, 1);
        end
        
        sesInfo.userInput.files = files;
        sesInfo.userInput.numOfSub = 1;
        sesInfo.userInput.numOfSess = length(aap.acq_details.selected_sessions);
        sesInfo.userInput.diffTimePoints = nscans;

        sesInfo.userInput.numOfPC1 = settings.PCA1;
        sesInfo.userInput.numOfPC2 = settings.PCA2;
        sesInfo.userInput.numComp = settings.numICs;
        sesInfo.userInput.preproc_type = settings.preproc;
        sesInfo.userInput.ICA_Options{10} = [settings.numICs];
        
        maskFile = aas_getfiles_bystream(aap, subj, 'firstlevel_brainmask');
        VM = spm_vol(maskFile);
        VM = VM(1);
        Y = spm_read_vols(VM);
        [path file ext] = fileparts(maskFile(1,:));
        VM.fname = fullfile(path, [file '.nii']);
        spm_write_vol(VM, Y); 
        sesInfo.userInput.maskFile = VM.fname;
        sesInfo = icatb_update_mask(sesInfo);

             
        save(fullfile(subjPath, parmFile) , 'sesInfo');
        aap = aas_desc_outputs(aap, subj, 'gift', fullfile(subjPath, parmFile) );
        
        
        icatb_runAnalysis(sesInfo, 1)
        
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end

end