function [aap resp]=aamod_emeg_forwardmodel(varargin)
% Create forward models
% Danny Mitchell 02/04/08

%% check task settings, subject, block etc
[aap subblock doit resp settings]=aa_emeg_checktasksettings(mfilename('fullpath'),varargin);
if ~doit; return; end

%% find files and decide whether to run task;
files=aas_emeg_findfiles(aap,settings.InputFilter,subblock);
if isempty(files); aas_log(aap,1,sprintf('\nFound no data! (Input filter is %s)\n',settings.InputFilter)); end

%% run for each file (except EEG or CDA)
for f=1:length(files);

    if ~isempty(regexp(files{f},'-eeg','ONCE')); continue; end
    if ~isempty(regexp(files{f},'_CDA_','ONCE')); continue; end
    if ~isempty(regexp(files{f},'_grp','ONCE')); continue; end
    if ~isempty(regexp(files{f},'_fused','ONCE')); continue; end

    %% load MEG header
    try rehash; load(files{f});
    catch
        fprintf('\nFailed to load file. Will try again in 10s in case there is an access conflict...\n');
        pause(10)
        try rehash; load(files{f},'-MAT');
        catch; debugnow
        end
    end

    %% check inv field and decide whether to overwrite
    try
        if isempty(D.inv) || settings.Overwrite==1; D.inv={struct([])}; end
    catch D.inv={struct([])};
    end;

    %% determine sensor type(s)
    fprintf('\nProcessing file %s',files{f});
    stype=regexp(files{f},'(-\w\w\w\w)\.','tokens');
    if isempty(stype); stype=''; % both sensor types
    else stype=stype{1}{1};
    end

    %% load subject meshes
    clear temp
    pth=fullfile(aap.acq_details.root,aap.acq_details.subjects(subblock{1}).megname,'structurals');
    try
        meshes=load(fullfile(pth,sprintf('Dinv%s.mat',stype)));
        meshes.D.inv{1}.mesh;
    catch
        try
            meshes=load(fullfile(pth,'Dinv.mat'));
            meshes.D.inv{1}.inverse.type='';
        catch error('Cortical mesh not found!');
        end
    end

    %% find the cell for this forward model if it exists else add new one
    tosave=false;
    D.val=1;
    for v=1:length(D.inv)
        try
            if strcmp(D.inv{v}.forward.method,settings.ForwardModel); D.val=v; break;
            elseif strcmp(D.inv{v}.forward.method,''); D.val=v; break; % type unspecified, so assume this inversion empty
            else D.val=length(D.inv)+1;
            end
        catch D.val=v; break % type unspecified, so assume this inversion empty
        end
    end

    if isempty(D.inv{D.val})
        D.inv{D.val}=meshes.D.inv{1};
        D.inv{D.val}.forward.method=settings.ForwardModel;
    end

    %% try to find cooregistration (this should be the same for any forward
    %% models in this file, so copy from first cell)
    try
        if isempty(D.inv{D.val}.datareg.eeg2mri) || settings.Overwrite==1;
            if D.val>1 && ~isempty(D.inv{1}.datareg.eeg2mri)
                D.inv{D.val}.datareg=D.inv{1}.datareg;
            end
        end
    catch % empty field might not have been created if using template MRI
        D.inv{D.val}.datareg.eeg2mri='';
    end

    %% Coregister MEG landmarks to MRI landmarks if necessary
    if isempty(D.inv{D.val}.datareg.eeg2mri)
        fprintf('\nCoregistering MEG with MRI...')
        D.inv{D.val}.datareg.sensors = D.channels.Loc';
        D.inv{D.val}.datareg.fid_eeg = D.channels.fid_eeg(:,1:3);
        D.inv{D.val}.datareg.fid_mri = D.inv{D.val}.mesh.smri_fids(:,1:3);
        D.inv{D.val}.datareg.headshape = D.channels.headshape(:,1:3);
        D.inv{D.val}.datareg.megorient = D.channels.Orient';
        D.inv{D.val}.datareg.scalpvert = D.inv{D.val}.mesh.tess_scalp.vert;
        D = spm_eeg_inv_datareg(D);
        checkdatareg(D,settings); % print diagnostic
        tosave=true;
    else fprintf('\nFound coregistration')
    end

    %% setup forward model (get leadfield matrix) if necessary
    D.inv{D.val}.method = 'Imaging';
    D.inv{D.val}.forward.surface2fit = 2; % inner skull
    gain_name = fullfile(D.path,[settings.ForwardModel stype '_SPMgainmatrix_' num2str(D.val) '.mat']);
    gxyz_name = strrep(gain_name,'matrix','matxyz');

    if ~exist(gain_name,'file');
        fprintf('\nConstructing forward model: ')
        addpath /imaging/local/Brainstorm/Toolbox % problem in /imaging/local/spm/spm5/meg_os.m at line 390; seems to derive from line 191.
        addpath /imaging/dm01/MEG/aaMEG % for overlapping_sphere.m graphics
        % my version allows to pass forward model type and stype as prefix for output filenames;
        % also show animated OS graphics. But let's go with the standard
        % version and rename files afterwards...
        %D = spm_eeg_inv_BSTfwdsol_dm(D,D.val,[settings.ForwardModel stype]);
        D = spm_eeg_inv_BSTfwdsol(D);
        eval(sprintf('!mv %s %s',D.inv{D.val}.forward.gainmat,gain_name))
        eval(sprintf('!mv %s %s',D.inv{D.val}.forward.gainxyz,gxyz_name))
        D.inv{D.val}.forward.gainmat=gain_name;
        D.inv{D.val}.forward.gainxyz=gxyz_name;
        tosave=true;
    else
        fprintf('\nFound forward model: %s', gain_name);
        if strcmp(D.inv{D.val}.forward.gainmat,gain_name) ...
        && strcmp(D.inv{D.val}.forward.gainxyz,gxyz_name);
            % aok
        else
            D.inv{D.val}.forward.gainmat=gain_name;
            D.inv{D.val}.forward.gainxyz=gxyz_name;
            tosave=true;
        end
    end

    if tosave; save(fullfile(D.path,D.fname),'D'); end

end % next file

return

%% plot and print diagnostic meshes of coregistration
function checkdatareg(D,settings)

outfile=fullfile(D.path,'figures','DataRegistration.jpg');
if exist(outfile,'file')==0 || settings.Overwrite==1
    warning off all; h = spm_figure('GetWin','Graphics'); spm_figure('Clear',h); warning on all;
    set(h,'renderer','opengl')
    fprintf('\nRendering images of coregistration...')

    try
        Lsens   = D.inv{1}.datareg.sens_coreg;
        Lfid    = D.inv{1}.datareg.fid_coreg;
        Lhsp    = D.inv{1}.datareg.hsp_coreg;
        Lfidmri = D.inv{1}.datareg.fid_mri;
    catch
        warndlg('please coregister these data')
        return
    end

    if size(Lsens,2)==6	% gradiometers, two coils
        for ch = 1:size(Lsens,1)
            if all(isfinite(Lsens(ch,4:6)))
                Lsens(ch,1:3) = (Lsens(ch,1:3)+Lsens(ch,4:6))/2;
            end
        end
        Lsens=unique(round(Lsens(:,1:3)),'rows');
    end

    for p=1:3

        axx=subplot(2,2,p);
        cla;

        if p==1 % render
            cortex=patch('vertices',D.inv{1}.mesh.tess_ctx.vert,'faces',D.inv{1}.mesh.tess_ctx.face,'EdgeColor','none','FaceColor',[0.8 1 1],'facealpha',1,'facelighting','phong');
            hold on
            skull=patch('vertices',D.inv{1}.mesh.tess_iskull.vert,'faces',D.inv{1}.mesh.tess_iskull.face,'EdgeColor','none','FaceColor','k','facealpha',0.2,'facelighting','gouraud');
            skin=patch('vertices',D.inv{1}.mesh.tess_scalp.vert,'faces',D.inv{1}.mesh.tess_scalp.face,'EdgeColor','none','FaceColor',[1 .7 .55],'facealpha',0.6,'facelighting','gouraud');

            material dull

            % now use Jason's code for 'cylinder' sensors:
            sth=2;  % Thickness of cylinders
            sr=5;   % Radius of cylinders
            sn=60;  % Number of elements in surface (large=pretty)
            H.sensors=[];
            for i=1:size(Lsens,1)
                x=Lsens(i,1); y=Lsens(i,2); z=Lsens(i,3);
                % This hack orients the sensors normal to (0,0,0):
                n=norm([x y z]);
                thx=acos(x/n); thy=acos(y/n); thz=acos(z/n);
                n1=n+sth;
                x1=n1*cos(thx); y1=n1*cos(thy); z1=n1*cos(thz);
                [H.sensors(i,1),H.sensors(i,2),H.sensors(i,3)]=Cylinder([x y z],[x1 y1 z1],sr,sn,'k',1,0);
            end
            set(H.sensors,'edgecolor','k','edgealpha',0.1,'facealpha',0.3);
            set(H.sensors(:,1),'facecolor',[.5 .5 .5]);     % cylinder body color
            set(H.sensors(:,2),'facecolor',[.5 .5 .5]);     % cylinder bottom surface color
            set(H.sensors(:,3),'facecolor',[.5 .5 .5]);     % cylinder top surface color
            set(H.sensors,...
                'AmbientStrength',         0.3,...
                'DiffuseStrength',         1,...
                'SpecularStrength',        0.7,...
                'SpecularExponent',        20,...
                'SpecularColorReflectance',1.0);

            % EEG fiducials or MEG coils (coreg.)
            h_fid   = plot3(Lfid(:,1),Lfid(:,2),Lfid(:,3),'o');
            set(h_fid,'MarkerFaceColor','g','MarkerSize',6,'MarkerEdgeColor','g');
            % MRI fiducials
            h_fidmr = plot3(Lfidmri(:,1),Lfidmri(:,2),Lfidmri(:,3),'o');
            set(h_fidmr,'MarkerFaceColor','r','MarkerSize',6,'MarkerEdgeColor','r');
            % headshape locations
            h_hsp   = plot3(Lhsp(:,1),Lhsp(:,2),Lhsp(:,3),'o');
            set(h_hsp,'MarkerFaceColor','b','MarkerSize',4,'MarkerEdgeColor','b');
            %         % Sensors (coreg.)
            %         h_sens  = plot3(Lsens(:,1),Lsens(:,2),Lsens(:,3),'o');
            %         set(h_sens,'MarkerFaceColor','none','MarkerSize',9,'MarkerEdgeColor',[0.4 0.4 0.4]);

            hold off
        else % just copy
            h_fid = copyobj(h_fid,axx);
            h_fidmr = copyobj(h_fidmr, axx);
            h_hsp = copyobj(h_hsp, axx);
            copyobj([H.sensors(:)', cortex, skull, skin],axx);
            axes(axx)
        end        
        axis image vis3d off;
        view([(p==2)*-90 (p==3)*90]); 
        camlight       
    end
    drawnow
    
    L=legend([h_fid(1) h_fidmr(1) h_hsp(1)], {'MEG fiducials','MRI fiducials','Digitised scalp points'}, ...
        'location','East');
    apos=get(axx,'position');
    lpos=get(L,'position');
    set(L,'position',[apos(1)+apos(3)+0.05 lpos(2:4)]);
    
    % convert to image
    f=getframe(h);
    clf; set(h,'renderer','painters');
    subplot('position',[0 0 1 1]);
    image(f.cdata); axis off image
    ann=annotation('textbox',[0 0 1 .06],'String',{outfile,sprintf('Date:...%s',datestr(now,0))},'edge','none');
    print('-djpeg90',outfile);
    delete(ann);
end
return