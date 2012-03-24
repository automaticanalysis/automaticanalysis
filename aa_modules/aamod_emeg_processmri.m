function [aap resp]=aamod_emeg_processmri(varargin)
% Copies MRI to subject's output directory, converting from dicom 
% if necessary (could let other aa functions do this?)
% Sets up D.inv structure for subject
% Calls spm_eeg_inv_spatnorm.m to normalise/segment (with initial bias correction)
% Calls spm_eeg_inv_meshing.m to normalise template meshes to subject head
% Automatically identifies fiducials on MRI if requested (not always well)
% Automatically masks out vitamin E capsules if requested
% Prints diagnostic images
%
% Danny Mitchell 13/03/08

%% check task settings, subject, block etc
[aap subblock doit resp settings]=aa_emeg_checktasksettings(mfilename('fullpath'),varargin);
if ~doit; return; end

%% get or create output directory and parameters structure

sub=aap.acq_details.subjects(subblock{1}).megname;
outpath=fullfile(aap.acq_details.root,sub,'structurals');
if exist(outpath,'dir')~=7; mkdir(outpath); end
outfile=fullfile(outpath,[sub '_T1.nii']);
    
params=fullfile(outpath,'Dinv.mat');
if ~exist(params,'file') || settings.Overwrite==1 
    load('/imaging/local/spm/spm5/EEGtemplates/defaults_eeg_inv.mat');
    D.inv{1} = invdefts.imag;
    D.inv{1}.mesh.sMRI = outfile;
else
    load('-mat',params);
end

D.inv{1}.mesh.canonical = 1;
D.inv{1}.mesh.Msize = settings.MeshSize; 
% spm meshes had [3004 4004 5004 7204] dipoles
% cbu_update meshes have [5124 8196 20484] dipoles

%% check whether structural exists
T1=aap.acq_details.subjects(subblock{1}).mriname;
if ~strcmp(T1,'missing')
    
    if ~exist(outfile,'file') && ~exist(strrep(outfile,'nii','img'),'file'); % don't overwrite: do this manually if needed
        
%% get filename of structural
        if ~exist(T1,'file') 
            %% look in central store
            fprintf('\nLooking for structural in central store...')
            filt=sprintf('^s%s.*$',aap.acq_details.subjects(subblock{1}).structuralfn);
            T1=spm_select('FPlist',aap.directory_conventions.centralstore_structurals,filt);
            T1=deblank(T1(1,:));
            if exist(T1,'file')~=2
                %% look for dicom directory
                fprintf('Failed\nLooking for directory of raw dicom files...')
                try
                    subpath=fullfile(aap.directory_conventions.rawdatadir,aap.acq_details.subjects(subblock{1}).mriname);
                    series=dir(subpath);
                    series=char(series.name);
                    ind=~cellfun('isempty',regexp(cellstr(series),'MPRAGE'));
                    T1=fullfile(subpath,deblank(series(find(ind,1),:)));
                    if exist(T1,'file')==7; fprintf('Found %s',T1);
                    else crash;
                    end
                catch
                   aas_log(aap,0,'\nFailed to find structural for %s\n', aap.acq_details.subjects(subblock{1}).megname);
                   return
                end
            else fprintf('Found %s',T1);
            end
        end
        
        if exist(T1,'file')==7; 
%% This is a directory: convert from dicom
            fprintf('\nLoading dicom headers...')
            dcms=spm_select('FPList',T1,'.*');
            hdrs=spm_dicom_headers(dcms);
            cd(outpath);
            fprintf('\nConverting to nifti...')
            spm_dicom_convert(hdrs,'all','flat','nii');
            dicomname=spm_select('List',pwd,'^.*nii$');
            dicomname=dicomname(cellfun('isempty',(regexp(cellstr(dicomname),'_T1.'))),:);
            movefile(dicomname,strrep(outfile,'.dcm','.nii'),'f');
        else
%% This is a file: copy to analysis directory
            fprintf('\nCopying to analysis directory...')
            cmd=sprintf('scp %s %s', T1, outfile);
            status=unix(cmd);
            if strcmp(T1(end-2:end),'img');
               % copy header too 
               cmd=sprintf('scp %s %s', strrep(T1,'img','hdr'), strrep(outfile,'img','hdr'));
               status=unix(cmd);
            end
        end;     
    end
    
%% mask out vitamin E capsules
    if settings.AutoMask
        if ~exist([outfile(1:end-4) '_mask.nii'],'file') || settings.Overwrite==1
            fprintf('\nMasking out vit. E capsules: ');
            aa_emeg_headmask(outfile,'display');
        end
        % Normalisation of masked image seems fine if structural is first
        % manually aligned. However with no manual realignment it seems
        % slightly less reliable than using the raw T1.
        % D.inv{1}.mesh.sMRI = [outfile(1:end-4) '_masked.nii'];
    end
    
%% Inverse normalise the canonical mesh to subject's native space
    [pth nam ext]=fileparts(outfile);
    biascorrected=fullfile(pth,['m' nam ext]);
    if ~strcmp(D.inv{1}.mesh.sMRI,biascorrected) || settings.Overwrite==1;
        %keyboard
        if ~exist(biascorrected,'file') || settings.Overwrite==1;
            % 1st pass normalisation for bias correction 
            fprintf('\nNormalising to template: First-pass bias correction...')
            estopts.regtype='';    % turn off affine:
            res= spm_preproc(D.inv{1}.mesh.sMRI,estopts);
            sn = spm_prep2sn(res);
            opts = struct('biascor',1,...
                  'GM',     [0 0 0],...
                  'WM',     [0 0 0],...
                  'CSF',    [0 0 0],...
                  'cleanup',0);
            spm_preproc_write(sn,opts);
        end
        D.inv{1}.mesh.sMRI=fullfile(pth,['m' nam ext]);
        save(params,'D');
    end
        
    if ~isfield(D.inv{1}.mesh,'wmMRI') || settings.Overwrite==1 
        % now redo the normalisation and invert the meshes 
        fprintf('\nSecond-pass of normalisation, and inversion of canonical meshes...\n')
        try D = spm_eeg_inv_meshing(D);	  
            % calls spm_eeg_inv_spatnorm
            % then spm_eeg_inv_getmasks
            % then loads and inverts canonical meshes 
        catch debugnow
        end
        
        % Replace the inv.norm.can. scalp and skull meshes with the subject's own
        % This should be more accurate as these can be segmanted easily (unlike cortex). 
        % But there could be a risk of the canonical cortex protruding through the
        % skull, although I haven't seen this happen yet...
        if settings.AutoMask
            D.inv{1}.mesh.msk_scalp=[outfile(1:end-4) '_mask.nii'];
        end
        D = spm_eeg_inv_getmeshes(D,2);
        save(params,'D');
    end 
    
%% automatically locate fiducials if requested
    if settings.AutoFindFiducials==1
        if ~isfield(D.inv{1}.mesh,'smri_fids') || settings.Overwrite==1
            if strcmpi('meg07_0081',sub) || strcmpi('meg07_0128',sub)
                D.lower='Correct for misdefined lateral fids';
            end
            D=FindFiducials(D);
            save(params,'D');
        end
    end
    
else %%%% this subject doesn't have a structural
    
%% generate template meshes without inverse normalising to subject
    if ~isfield(D.inv{1}.mesh.tess_iskull,'face') || settings.Overwrite==1 
        fprintf('No MRI found; using template mesh\n')
        D = spm_eeg_inv_template(D);
        save(params,'D');
    end
    
%% automatically define fiducials if requested
    if settings.AutoFindFiducials==1
        if strcmpi('meg07_0081',sub) || strcmpi('meg07_0128',sub)
            D.inv{1}.mesh.wmri_fids=[0,81,-50;-85,-21,-70;80,-19,-70];
        else
            D.inv{1}.mesh.wmri_fids=[0,81,-50;-85,-25,-53;80,-22,-53];
        end
        D.inv{1}.mesh.smri_fids=D.inv{1}.mesh.wmri_fids;
        save(params,'D');
    end
    
end

%% prepare figure
plotfile=spm_select('Cpath',fullfile(outpath,'..','figures',[sub '_meshes_MRI.jpg']));
if ~exist(plotfile,'file') || settings.Overwrite==1
    h = spm_figure('GetWin','Graphics'); spm_figure('Clear',h)
    ann=annotation('textbox',[0.6 .94 0.4 .06],'String', ...
        {sprintf('Subject: %s',aap.acq_details.subjects(subblock{1}).megname), ...
        sprintf('Date:...%s',datestr(now,0))}, ...
        'edge','none');
    
%% show estimated fiducial locations in subject's MRI space:
    global st
    offmm=D.inv{D.val}.mesh.smri_fids;
    brain=D.inv{1}.mesh.sMRI;
    spm_orthviews('reset'); 
    spm_orthviews('image',brain,[0 0.01 0.5 0.33]);
    spm_orthviews('reposition',offmm(1,:)); st.vols{1}=[];xlabel('Nasion')
    spm_orthviews('image',brain,[0 0.35 0.5 0.33]);
    spm_orthviews('reposition',offmm(2,:)); st.vols{1}=[];xlabel('LeftEar')
    spm_orthviews('image',brain,[0 0.69 0.5 0.33]);
    spm_orthviews('reposition',offmm(3,:)); st.vols{1}=[];xlabel('RightEar')

%% show meshes and normalised structural
    set(h,'renderer','opengl')
    for p=1:4
        subplot(4,2,p*2)
        if p==1 % normalisation
            axis image off
            try
                V=spm_vol(D.inv{1}.mesh.wmMRI);
                Y=spm_read_vols(V);
                imagesc(rot90(squeeze(Y(floor(end/2),:,:))));
                title('Normalised structural')
                clear Y V
                axis equal off
            catch
                title('No MRI. Using template mesh.')
            end
        else % meshes
            cla;
            cortex=patch('vertices',D.inv{1}.mesh.tess_ctx.vert,'faces',D.inv{1}.mesh.tess_ctx.face,'EdgeColor','none','FaceColor',[0.8 1 1],'facelighting','phong');
            hold on
            skull=patch('vertices',D.inv{1}.mesh.tess_iskull.vert,'faces',D.inv{1}.mesh.tess_iskull.face,'EdgeColor','none','FaceColor','k','facealpha',0.4,'facelighting','gouraud');
            skin=patch('vertices',D.inv{1}.mesh.tess_scalp.vert,'faces',D.inv{1}.mesh.tess_scalp.face,'EdgeColor','none','FaceColor',[1 .7 .55],'facealpha',0.5,'facelighting','gouraud');
            axis image vis3d off;
            view((p==2)*90,(p==3)*90);
            camlight
            material dull
            hold off
            drawnow
        end
    end

    fprintf('Printing...'); print('-djpeg90',plotfile);
    set(h,'renderer','painters')
    delete(ann)
end

return