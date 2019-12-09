% AA module
% Runs BET (FSL Brain Extration Toolbox) on structural (usually)
% [For best functionality, it is recommended you run this after
% realignment and before writing the normalised EPI image
% If you do it before estimating the normalisation, make sure you normalise
% to a scull-stripped template, if at all possible!]

function [aap,resp]=aamod_bet(aap,task,varargin)

resp='';

switch task
    case 'report'
        subj = varargin{1};
        domain = aap.tasklist.currenttask.domain;
        localpath = aas_getpath_bydomain(aap,domain,cell2mat(varargin));

        if isempty(spm_select('List',localpath,'diagnostic_.*.jpg'))
            diag(aap,varargin{:});
        end
        
        for fdiag = cellstr(spm_select('FPList',localpath,'diagnostic_.*.jpg'))'
            aap = aas_report_add(aap,subj,'<table><tr><td>');
            aap=aas_report_addimage(aap,subj,fdiag{1});
            [p, f] = fileparts(fdiag{1}); avipath = fullfile(p,[strrep(f(1:end-2),'slices','avi') '.avi']);
            if exist(avipath,'file'), aap=aas_report_addimage(aap,subj,avipath); end
            aap = aas_report_add(aap,subj,'</td></tr></table>');
        end
    case 'doit'
        % Whatever domain we're operating at must have subj as top domain
        fslext=aas_getfslext(aap);
        
        % Get inputs
        streamname=aap.tasklist.currenttask.inputstreams.stream{1};
        % Let us use the native space...
        Simg = aas_getfiles_bystream(aap,aap.tasklist.currenttask.domain,[varargin{:}],streamname);
        
        % Which file is considered, as determined by the structural parameter!
        if size(Simg,1) > 1
            Simg = deblank(Simg(aap.tasklist.currenttask.settings.structural, :));
            aas_log(aap,false,sprintf('\tWARNING: Several structurals found, considering: %s', Simg))
        end
        
        [pth nme ext]=fileparts(Simg);
        
        outStruct=fullfile(pth,['bet_' nme fslext]);

        fslcommand = sprintf('bet %s %s -f %f -v',Simg,outStruct, ...
                aap.tasklist.currenttask.settings.bet_f_parameter);

        if aap.tasklist.currenttask.settings.robust
            aas_log(aap,false,'1st BET pass (recursive) to find optimal centre of gravity and radius')
            fslcommand = [fslcommand ' -R'];
        else aap.tasklist.currenttask.settings.biasneck
            aas_log(aap,false,'1st BET pass (using bias field and neck cleanup) to find optimal centre of gravity and radius')
            fslcommand = [fslcommand ' -B'];
        end
        
        [junk, w]=aas_runfslcommand(aap, fslcommand);
        
        aas_log(aap,false,sprintf('Bet output: %s',w));
        
        % This outputs last radius from recursive command...
        indxS = strfind(w, 'radius');
        indxS = indxS(end) + 7;
        indxE = strfind(w(indxS:end), ' mm');
        indxE = indxE(1) - 2;
        SRad = w(indxS:indxS+indxE);
        
        % We don't extract the centre of gravity from here, since it needs
        % to be input in voxels... Instead get it from betted image
        [V Y] = aas_spm_vol(outStruct);
        Y = Y > 0;
        indY = find(Y);
        [subY_x subY_y subY_z] = ind2sub(size(Y), indY);
        COG = [mean(subY_x), mean(subY_y), mean(subY_z)];
        
        aas_log(aap,false,sprintf('\t...calculated c-o-g (vox): %0.4f %0.4f %0.4f  and radius (mm): %s', ...
            COG(1), COG(2), COG(3), SRad))
        
        if aap.tasklist.currenttask.settings.masks
            aas_log(aap,false,'2nd BET pass extracting brain masks')
            % Run BET [-A Now let's get the brain masks and meshes!!]
            [junk, w]=aas_runfslcommand(aap, ...
                sprintf('bet %s %s -f %f -c %0.4f %0.4f %0.4f -r %s -v -A',Simg,outStruct, ...
                aap.tasklist.currenttask.settings.bet_f_parameter, COG(1), COG(2), COG(3), SRad)...
                );
            
            % This outputs last radius from recursive command...
            indxS = strfind(w, 'radius');
            indxS = indxS(end) + 7;
            indxE = strfind(w(indxS:end), ' mm');
            indxE = indxE(1) - 2;
            SRad = w(indxS:indxS+indxE);
            
            % We don't extract the centre of gravity from here, since it needs
            % to be input in voxels... Instead get it from betted image
            [V Y] = aas_spm_vol(outStruct);
            Y = Y > 0;
            indY = find(Y);
            [subY_x subY_y subY_z] = ind2sub(size(Y), indY);
            COG = [mean(subY_x), mean(subY_y), mean(subY_z)];
            
            aas_log(aap,false,sprintf('\t...final c-o-g (vox): %0.4f %0.4f %0.4f  and radius (mm): %s', ...
                COG(1), COG(2), COG(3), SRad))
        end
        
        % BRAIN MASK (slightly different from inskull)
        
        [V mY] = aas_spm_vol(fullfile(pth,['bet_' nme fslext]));
        
        % Mask out non-brain
        mY = mY > 0;
        
        % Then write out actual BET mask
        V.fname = fullfile(pth, ['bet_' nme '_brain_mask' ext]);
        spm_write_vol(V,mY);
        
        % Get the mask images
        D = dir(fullfile(pth, 'bet*mask*'));
        outMask = '';
        for d = 1:length(D)
            outMask = strvcat(outMask, fullfile(pth, D(d).name));
        end
        
        % Get also the meshes
        D = dir(fullfile(pth, 'bet*mesh*'));
        outMesh = '';
        for d = 1:length(D)
            outMesh = strvcat(outMesh, fullfile(pth, D(d).name));
        end
        
        %% DESCRIBE OUTPUTS!
        
        % Structural image after BETting
        outstreams = aas_getstreams(aap,'output'); % image, mask, mesh
        aap=aas_desc_outputs(aap,aap.tasklist.currenttask.domain,[varargin{:}],streamname,outStruct);
        aap=aas_desc_outputs(aap,aap.tasklist.currenttask.domain,[varargin{:}],outstreams{2},outMask);
        if aap.tasklist.currenttask.settings.masks
            aap=aas_desc_outputs(aap,aap.tasklist.currenttask.domain,[varargin{:}],outstreams{3},outMesh);
        end
        
        try diag(aap,varargin{:}); catch, end
end
end

function diag(aap,varargin)
%% DIAGNOSTIC IMAGE
%% Draw structural image...
inpstreams = aas_getstreams(aap,'input'); % 
outstreams = aas_getstreams(aap,'output'); % 
Simg = aas_getfiles_bystream(aap,aap.tasklist.currenttask.domain,[varargin{:}],inpstreams{1},'input');
outStruct = aas_getfiles_bystream(aap,aap.tasklist.currenttask.domain,[varargin{:}],outstreams{1},'output');
outMesh = aas_getfiles_bystream(aap,aap.tasklist.currenttask.domain,[varargin{:}],outstreams{3},'output');

spm_check_registration(Simg)

% This will only work for 1-7 masks
OVERcolours = {[1 0 0], [0 1 0], [0 0 1], ...
    [1 1 0], [1 0 1], [0 1 1], [1 1 1]};

indx = 0;

% Colour the brain extracted bit pink
spm_orthviews('addcolouredimage',1,outStruct, [0.9 0.4 0.4])
% Add mesh outlines, to see if BET has worked properly!
if aap.tasklist.currenttask.settings.masks
    for r = 1:size(outMesh,1)
        if strfind(outMesh(r,:), aas_getfslext(aap))
            indx = indx + 1;
            spm_orthviews('addcolouredimage',1,outMesh(r,:), OVERcolours{indx})
        end
    end
end
%% Diagnostic VIDEO of masks
aas_checkreg_avi(aap, {aap.tasklist.currenttask.domain cell2mat(varargin)}, 0)

spm_orthviews('reposition', [0 0 0])

try figure(spm_figure('FindWin', 'Graphics')); catch; figure(1); end;
set(gcf,'PaperPositionMode','auto','Renderer','zbuffer');
print('-djpeg','-r75',fullfile(aas_getpath_bydomain(aap,aap.tasklist.currenttask.domain,[varargin{:}]), ...
    ['diagnostic_' aap.tasklist.main.module(aap.tasklist.currenttask.modulenumber).name '_' aap.acq_details.subjects(varargin{1}).subjname '.jpg']));
end
