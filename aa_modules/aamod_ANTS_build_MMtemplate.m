% AA module - build template from a set of structural scans, using ANTS
% It will build a multimodal template, based on the various streams that it
% receives... (e.g. structurals, segmentations, rois, etc.)

function [aap,resp]=aamod_ANTS_build_MMtemplate(aap,task)

resp = '';

switch task
    case 'doit'
        
        % Make template path
        Tpth = fullfile(aas_getstudypath(aap), 'ANTStemplate');
        if ~exist(Tpth, 'dir')
            mkdir(Tpth)
        else
            unix(['rm -rf ' Tpth])
            mkdir(Tpth)
        end
        
        % Expects only one input stream...
        streams=aap.tasklist.currenttask.inputstreams.stream;
        
        fid = fopen(fullfile(Tpth, 'MMtemplate.txt'), 'w');
        for subj = 1:length(aap.acq_details.subjects)
            % We must make a .txt file, with each line containing all the
            % images for a subject...
            if subj ~= 1
                fprintf(fid, '\n');
            end
            
            if subj == 1
                modalities = 0;
            end
            
            for st = 1:length(streams)
                
                Simg = aas_getfiles_bystream(aap,subj,streams{st});
                
                if isfield(aap.tasklist.currenttask.settings, [streams{st} '_images'])
                    Simg = Simg(aap.tasklist.currenttask.settings.([streams{st} '_images']), :);
                end
                
                for ts = 1:size(Simg,1)
                    if subj == 1
                        modalities = modalities + 1;
                    end
                    
                    % Fileparts to get extension of file...
                    [junk, junk, Sext] = fileparts(deblank(Simg(ts,:)));
                    
                    fn = fullfile(Tpth, sprintf('subj%04d_%d_%d%s', subj, st, ts, Sext));
                    
                    % Write to .txt...
                    fprintf(fid, [fn ' ']);
                    % Copy images to right location
                    copyfile(Simg(ts,:), fn);
                end
            end
        end
        
        %% Use ANTS to make a template!
        
        % Set the ANTS path
        setenv('ANTSPATH', fullfile(aap.directory_conventions.ANTSdir, 'bin/'))
        ANTSpath = [' sh ' fullfile(getenv('ANTSPATH'), 'antsMultivariateTemplateConstruction.sh') ' '];
        % Add the path with functions to interact with torque (qsub) <-- changed
        %setenv('PATH', [getenv('PATH') ':' fullfile(aap.directory_conventions.fieldtripdir, 'qsub')])
        
        % What we get out...
        outfiles = '-o ANTS ';
        
        % Dimension number (always 3 for structural)
        Ndim = '-d 3 ';
        
        options = aap.tasklist.currenttask.settings.extraoptions;
        
        ANTS_command = [ ANTSpath Ndim outfiles options ' -k ' num2str(modalities) ' MMtemplate.txt'];
        
        cd(Tpth)
        
        % Run ANTS
        fprintf('Running ANTS using command:\n')
        fprintf([ANTS_command '\n'])
        
        [s w] = aas_shell(ANTS_command);
        disp(w)
        
        %% Describe the outputs
        outTemp = '';
        for t = 1:modalities
            unix(['gunzip ' fullfile(Tpth, ['ANTStemplate' num2str(t-1) '.nii.gz'])])
            outTemp = strvcat(outTemp, fullfile(Tpth, ['ANTStemplate' num2str(t-1) '.nii']));
        end
        aap = aas_desc_outputs(aap,'ANTStemplate', outTemp);
        
        % Delete other things
        delete(fullfile(Tpth,'*nii.gz'))
        delete(fullfile(Tpth,'*txt'))
        delete(fullfile(Tpth,'subj*'))
        delete(fullfile(Tpth,'rigid*'))
        D = dir(Tpth);
        for d = 3:length(D)
            if isdir(fullfile(Tpth, D(d).name))
                rmdir(fullfile(Tpth, D(d).name), 's')
            end
        end
        
        %% Diagnostic image?
        mriname = aas_prepare_diagnostic(aap,subj);
                
        %% Draw templates
        
        spm_check_registration(outTemp)
        
        spm_orthviews('reposition', [0 0 0])
        
        try figure(spm_figure('FindWin', 'Graphics')); catch; figure(1); end;
        print('-djpeg','-r150',fullfile(aap.acq_details.root, 'diagnostics', ...
            [mfilename '__' mriname '.jpeg']));
end