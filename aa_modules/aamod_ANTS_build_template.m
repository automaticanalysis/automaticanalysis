% AA module - build template from a set of structural scans, using ANTS

function [aap,resp]=aamod_ANTS_build_template(aap,task)

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
        streams=aap.tasklist.currenttask.inputstreams;
        
        for subj = 1:length(aap.acq_details.subjects)
            Simg = aas_getfiles_bystream(aap,subj,streams.stream{:});
            
            if size(Simg,1) > 1
                Simg = Simg(aap.tasklist.currenttask.settings.structural,:);
                aas_log(aap, false, sprintf('Found more than 1 image, using %s', ...
                    Simg));
            end
            
            % Fileparts to get extension of file...
            [junk, junk, Sext] = fileparts(Simg);
            
            % Copy images to right location
            copyfile(Simg, fullfile(Tpth, sprintf('subj%04d%s', subj, Sext)));
        end
        
        %% Use ANTS to make a template!
       
        % Set the ANTS path
        setenv('ANTSPATH', fullfile(aap.directory_conventions.ANTSdir, 'bin/'))
        ANTSpath = [' sh ' fullfile(getenv('ANTSPATH'), 'buildtemplateparallel.sh') ' '];
        % Add the path with functions to interact with torque (qsub) <-- changed
        %setenv('PATH', [getenv('PATH') ':' fullfile(aap.directory_conventions.fieldtripdir, 'qsub')])
        
        % What we get out...
        outfiles = '-o ANTS ';
        
        % Dimension number (always 3 for structural)
        Ndim = '-d 3 ';
        
        options = aap.tasklist.currenttask.settings.extraoptions;
        
        ANTS_command = [ ANTSpath Ndim outfiles options ' subj*'];
        
        cd(Tpth)
        
        % Run ANTS
        fprintf('Running ANTS using command:\n')
        fprintf([ANTS_command '\n'])
        
        [s, w] = aas_shell(ANTS_command);
        
        %% Describe the outputs
        unix(['gunzip ' fullfile(Tpth, 'ANTStemplate.nii.gz')])
        outTemp = fullfile(Tpth, 'ANTStemplate.nii');
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
        
        %% DIAGNOSTIC
        subjname = aas_prepare_diagnostic(aap, subj);
        
        % Draw template
        
        spm_check_registration(outTemp)
        
        spm_orthviews('reposition', [0 0 0])
        
        f = spm_figure('GetWin','Graphics');
        set(f,'Renderer','zbuffer');
        print(f, '-djpeg','-r150',fullfile(aap.acq_details.root, 'diagnostics', ...
            [mfilename '__' subjname '.jpeg']));
end