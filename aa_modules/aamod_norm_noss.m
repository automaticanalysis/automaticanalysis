% AA module - normalisation using normalise or two pass procedure with segment
% [aap,resp]=aamod_norm_noss(aap,task,subj)
% Depending on aap.tasksettings.aamod_norm_noss.usesegmentnotnormalise
% If 0 - classic Normalisation
% If 1 - use Segment with two pass procedure: a first pass to correct
%  for intensity variation across our structural images (probably due to
%  inhomogeneous SNR across space of 12 channel coil); and a second pass
%  to then do the segmentation
% _noss version does not use skull stripping
% subj=subject num
% Rhodri Cusack & Daniel Mitchell MRC CBU 2006
% based on originals by Rik Henson, Matthew Brett

function [aap,resp]=aamod_norm_noss(aap,task,subj,sess)
resp='';

switch task
    case 'report' % [TA]
        domain=aap.tasklist.currenttask.domain;
        switch domain
            case 'subject'
                indices=subj;
            case 'diffusion_session'
                indices=[subj sess];
        end;
        if isempty(dir(fullfile(aas_getsubjpath(aap,subj),'diagnostic_*.jpg')))
            diag(aap,domain,indices);
        end
        if strcmp(domain,'subject')
            fdiag = dir(fullfile(aas_getsubjpath(aap,subj),'diagnostic_*.jpg'));
            for d = 1:numel(fdiag)
                aap = aas_report_add(aap,subj,'<table><tr><td>');
                imgpath = fullfile(aas_getsubjpath(aap,subj),fdiag(d).name);
                aap=aas_report_addimage(aap,subj,imgpath);
                [p f] = fileparts(imgpath); avipath = fullfile(p,[strrep(f(1:end-2),'slices','avi') '.avi']);
                if exist(avipath,'file'), aap=aas_report_addimage(aap,subj,avipath); end
                aap = aas_report_add(aap,subj,'</td></tr></table>');
            end
        end;
    case 'doit'
        
        
        settings = aap.tasklist.currenttask.settings;
        
        defs =aap.spm.defaults.normalise;
        defs.estimate.weight = '';
       
        % Template image; here template image not skull stripped
        % [AVG] Changed to allow specification of any T1 template, does not
        % need to be in the SPM folder any more...
        temp_imgs = aap.directory_conventions.T1template;
        if (~exist(temp_imgs,'file')) % try in SPM
            if temp_imgs(1) ~= '/', temp_imgs = fullfile(fileparts(which('spm')),temp_imgs); end
        end
        if ~exist(temp_imgs, 'file')
            aas_log(aap, true, sprintf('Couldn''t find template T1 image %s.', temp_imgs));
        end
        
        clear imgs;
        
        %% Get structural
        % [AVG] Modified the way we get the structural, to be more aa4-like
        inStream = aap.tasklist.currenttask.inputstreams(1).stream{1};

        domain=aap.tasklist.currenttask.domain;
        
        if strcmp(domain,'diffusion_session')
            indices=[subj sess];
            Simg = aas_getfiles_bystream(aap,'diffusion_session',[subj sess],inStream);
        end;
        
        try 
            Simg = aas_getfiles_bystream(aap,subj,inStream);
            domain='subject';
            indices=[subj];
        catch
            Simg = aas_getfiles_bystream(aap,subj,1,inStream);
            domain='session';
            indices=[subj 1];
        end
        
        % Cheap and cheerful way of ensuring only one file is considered!
        if size(Simg,1) > 1
            Simg = deblank(Simg(1,:));
            aas_log(aap,0,sprintf('Found more than one structural so using:\n%s',Simg));
        end
        % Get structural directory for this subject
        [Spth, Sfn, Sext] = fileparts(Simg);
        
        %% Set up normalisation, etc.
        
        % Set the mask for subject to empty by default
        objMask = ''; % object mask
        
        % Check if there is an affine starting estimate for this subject
        allSubj = strcmp({settings.subject(:).name}, '*');
        thisSubj = strcmp({settings.subject(:).name}, aap.acq_details.subjects(subj).subjname);
        
        if any(allSubj)
            subjectOptions = settings.subject(allSubj);
        elseif any(thisSubj)
            subjectOptions = settings.subject(thisSubj);
        else
            subjectOptions = settings.subject(1);
        end   
        
        % Because we are going reslice later (with undistort_reslice)
        % We don't reslice anything except the image to be normalized
        
        % call the SPM segment or normalize function to do the work
        if (aap.tasklist.currenttask.settings.usesegmentnotnormalise)
            % 2 stage process, as proposed by RH, to increased robustness [djm 13/03/06]
            
            %%%%%%%% 1st pass:
            aas_log(aap,false,'Running first pass of norm_noss (get bias corrected structural)')
            estopts.regtype='';    % turn off affine:
            if ~isfield(estopts,'fudge'), estopts.fudge = 5; end % compatiblity
            out = spm_preproc(Simg, estopts);
            [sn,isn]   = spm_prep2sn(out);
            
            % only write out attenuation corrected image
            writeopts.biascor = 1;
            writeopts.GM  = [0 0 0];
            writeopts.WM  = [0 0 0];
            writeopts.CSF = [0 0 0];
            writeopts.cleanup = [0];
            spm_preproc_write(sn, writeopts);
            
            %%%%%%%% 2nd pass using attenuation corrected image
            aas_log(aap,false,'Running second pass of norm_noss (estimate normalisation)')
            % mstruc should be the attenuation corrected image
            % look for m prefixed filename
            mSimg = fullfile(Spth,['m' Sfn Sext]);
            if size(mSimg,1)>1
                aas_log(aap,0,sprintf('Found more than one attenuated structural so using first:\n%s',Simg));
            end
            
%             estopts.regtype='mni';    % turn on affine again
            
            % [TA] SPM12 compatibility
            try estopts = aap.spm.defaults.preproc; catch, estopts = aap.spm.defaults.old.preproc; end
            if iscell(estopts.tpm)
                estopts.tpm = char(estopts.tpm);
            end
            
            if isfield(aap.tasklist.currenttask.settings,'estopts')
                fields = fieldnames(aap.tasklist.currenttask.settings.estopts);
                for f = 1:numel(fields)
                    if ~isempty(aap.tasklist.currenttask.settings.estopts.(fields{f}))
                        estopts.(fields{f}) = aap.tasklist.currenttask.settings.estopts.(fields{f});
                    end
                end
            end
            
            % Load header of image to be normalized
            V=spm_vol(mSimg);
            
            % Now adjust parameters according to starting offset parameters
            % as requested in task settings
            startingParameters=[0 0 0   0 0 0   1 1 1   0 0 0];
            ParameterFields={'x','y','z', 'pitch','roll','yaw', 'xscale','yscale','zscale', 'xaffign','yaffign','zaffign'};

            if ~isempty(subjectOptions.affineStartingEstimate)
                aSE = subjectOptions.affineStartingEstimate;
            elseif ~isempty(aap.tasklist.currenttask.settings.affinestartingestimate)
                aSE = aap.tasklist.currenttask.settings.affinestartingestimate;
                aas_log(aap, 0, 'Warning: the option <affinestartinestimate> will soon be removed from aamod_norm_noss. Please use the <subject> option to apply affine starting estimates');
            else
                aSE = [];
            end
            
            if isstruct(aSE)
                fnames = fieldnames(aSE);
                for fieldind=1:length(fnames)
                    % Which element in StartingParameters does this refer to?
                    whichitem=find([strcmp(fnames{fieldind},ParameterFields)]);
                    % Generate a helpful error if it isn't recognised
                    if (isempty(whichitem))
                        err=sprintf('Unexpected field %s in header file aamod_norm_noss.xml - expected one of ',fnames{fieldind});
                        err=[err sprintf('%s\t',ParameterFields)];
                        aas_log(aap,true,err);
                    end
                    % Put this in its place
                    if ~isempty(aSE.(fnames{fieldind}))
                        startingParameters(whichitem)=aSE.(fnames{fieldind});
                    end;
                end;
            else
               startingParameters = aSE; 
            end

            %[AVG] Save original V.mat parameters
            oldMAT = V.mat;
            
            % Adjust starting orientation of object image as requested
            startingAffine=spm_matrix(startingParameters);
            V.mat=startingAffine*V.mat;
            
            % Run normalization
            out = spm_preproc(V,estopts);
            
            % Adjust output Affine to reflect fiddling of starting
            % orientation of object image
            out.Affine=out.Affine*startingAffine;
            
            % [AVG] Instead we set the out.image parameters to our original
            % structural image!
            out.image.mat = oldMAT;
            
            [sn,isn]   = spm_prep2sn(out);
            
            % [AVG] DEBUG:
            % We could print out the Affines in an orderly way, so that
            % the experimenter can see what is being printed out...
            %{
            matstr = strrep_multi(mat2str(oldMAT,2),{';', ' '},{'\n' '\t'}); matstr = matstr(2:end-1);
            aas_log(aap,false,sprintf('Initial structural Affine\n%s',matstr))
            matstr = strrep_multi(mat2str(StartingAffine,2),{';', ' '},{'\n' '\t'}); matstr = matstr(2:end-1);
            aas_log(aap,false,sprintf('Initial Affine transform\n%s',matstr))
            matstr = strrep_multi(mat2str(out.Affine,2),{';', ' '},{'\n' '\t'}); matstr = matstr(2:end-1);
            aas_log(aap,false,sprintf('Out image Affine\n%s',matstr))
            matstr = strrep_multi(mat2str(sn.Affine,2),{';', ' '},{'\n' '\t'}); matstr = matstr(2:end-1);
            aas_log(aap,false,sprintf('Spatial Normalisation Affine\n%s',matstr))
            %}
            
            aas_log(aap,false,'Writing out the segmented images')
            % write out GM , WM, CSF native + unmod normalised
            writeopts.biascor = 1;
            writeopts.GM  = [0 1 1];
            writeopts.WM  = [0 1 1];
            writeopts.CSF = [0 1 1];
            writeopts.cleanup = [0];
            spm_preproc_write(sn,writeopts);
            
            SNmat = fullfile(Spth, [Sfn '_seg_sn.mat']);
            invSNmat = fullfile(Spth, [Sfn '_seg_inv_sn.mat']);
            savefields(SNmat,sn);
            savefields(invSNmat,isn);
            aap=aas_desc_outputs(aap,subj,'normalisation_seg_sn',SNmat);
            aap=aas_desc_outputs(aap,subj,'normalisation_seg_inv_sn',invSNmat);
            
            % [AVG] this output is completely different from .xml
            % [RC] I like having separate streams, as it removes the need
            % for later file filtering, which is package dependent
            % So lets support both...
            tiss={'grey','white','csf'};
            [mSpth mSnme mSext]=aas_fileparts(mSimg);
            for tissind=1:3
                aap=aas_desc_outputs(aap,subj,sprintf('tissue_%s',tiss{tissind}),fullfile(Spth,sprintf('wc%d%s',tissind,[mSnme mSext])));
                aap=aas_desc_outputs(aap,subj,sprintf('unmod_tissue_%s',tiss{tissind}),fullfile(Spth,sprintf('c%d%s',tissind,[mSnme mSext])));
            end
            
            % [AVG] also group it all into segmentation stream
            outSeg = '';
            d = 0;
            while ~isnan(d)
                d = d+1;
                if exist(fullfile(Spth,sprintf('c%d%s',d,['m' Sfn Sext])), 'file')
                    outSeg = strvcat(outSeg, fullfile(Spth,sprintf('c%d%s',d,['m' Sfn Sext])));
                    outSeg = strvcat(outSeg, fullfile(Spth,sprintf('wc%d%s',d,['m' Sfn Sext])));
                else
                    d = NaN;
                end
            end
            aap=aas_desc_outputs(aap,domain,indices,'segmentation',outSeg);
            
            % [TA] replace the structural with the bias-corrected one
            Simg = fullfile(Spth,['mm' Sfn Sext]);
            Sout = fullfile(Spth,[aap.spm.defaults.normalise.write.prefix 'mm' Sfn Sext]);
        else
            % Make the default normalization parameters file name
            % Turn off template weighting
            % SPM defaults
            SNmat = fullfile(Spth, [Sfn '_sn.mat']);
            spm_normalise(temp_imgs, Simg, SNmat,...
                defs.estimate.weight, objMask, ...
                defs.estimate);
            aap=aas_desc_outputs(aap,subj,'normalisation_seg_sn',SNmat);
            % SPM2 normalization doesn't generate the inverse transformation
            %             invSNmat = fullfile(Spth, [Sfn '_seg_inv_sn.mat']);
            % aap=aas_desc_outputs(aap,subj,'normalisation_seg_inv_sn',invSNmat);
            Sout = fullfile(Spth,[aap.spm.defaults.normalise.write.prefix Sfn Sext]);
        end
        
        spm_write_sn(Simg,SNmat,defs.write);
        % [AVG] we need to add all the outputs, including warped structural
        % [AVG] It is probably best to save the 2nd bias-corrected image
        
        % [CW] But we don't have a bias corrected image if we didn't use
        % segmentation.

        % convert structural (bias-corrected) image into int16 (for FSL diag)
        V = spm_vol(Sout); Y = spm_read_vols(V);
        V.dt = [spm_type('int16') spm_platform('bigend')];
        spm_write_vol(V,Y);
        
        aap=aas_desc_outputs(aap,domain,indices,'structural', strvcat(Simg, Sout));

        if settings.diagnostic && strcmp(aap.options.wheretoprocess,'localsingle')
            diag(aap,domain,indices);
        end
        
    case 'checkrequirements'
        % Template image; here template image not skull stripped
        T1file = aap.directory_conventions.T1template;
        if (~exist(T1file,'file')) % try in SPM
            if T1file(1) ~= '/', T1file = fullfile(fileparts(which('spm')),T1file); end
        end
        if (~exist(T1file,'file'))
            aas_log(aap,true,sprintf('T1 template file %s not found, check aap.directory_conventions.T1template\n',T1file));
        end
end
end
%------------------------------------------------------------------------
function savefields(fnam,p)
if length(p)>1, aas_log([],true,'Can''t save fields.'); end
if str2double(version('-release'))>=14,
    save(fnam,'-V6','-struct','p');
else
    save(fnam,'-struct','p');
end
end
%------------------------------------------------------------------------
function diag(aap,domain,indices) % [TA]
% SPM, AA
Simg = aas_getfiles_bystream(aap,domain,indices,'structural','output'); Simg = Simg(1,:);
localpath = aas_getpath_bydomain(aap,domain,indices);
outSegAll = aas_getfiles_bystream(aap,domain,indices,'segmentation','output');

outSeg = outSegAll([1 3 5],:);
outNSeg = outSegAll([2 4 6],:);

OVERcolours = {[1 0 0], [0 1 0], [0 0 1]};

%% Draw native template
spm_check_registration(Simg)
% Add normalised segmentations...
for r = 1:size(outSeg,1)
    spm_orthviews('addcolouredimage',1,deblank(outSeg(r,:)), OVERcolours{r});
end

spm_orthviews('reposition', [0 0 0])

try f = spm_figure('FindWin', 'Graphics'); catch; f = figure(1); end;
set(f,'Renderer','zbuffer');
print(f,'-djpeg','-r150',...
    fullfile(localpath,['diagnostic_' aap.tasklist.main.module(aap.tasklist.currenttask.modulenumber).name '_N_blob.jpg']));

%% Draw warped template
tmpfile = fullfile(aap.directory_conventions.fsldir,'data','standard','MNI152_T1_1mm.nii.gz'); % use FSL highres
gunzip(tmpfile,localpath); tmpfile = fullfile(localpath,'MNI152_T1_1mm.nii');

spm_check_registration(tmpfile)
% Add normalised segmentations...
for r = 1:size(outNSeg,1)
    spm_orthviews('addcolouredimage',1,outNSeg(r,:), OVERcolours{r})
end
spm_orthviews('reposition', [0 0 0])

try f = spm_figure('FindWin', 'Graphics'); catch; f = figure(1); end;
set(f,'Renderer','zbuffer');
print(f,'-djpeg','-r150',...
    fullfile(localpath,['diagnostic_' aap.tasklist.main.module(aap.tasklist.currenttask.modulenumber).name '_W_blob.jpg']));

close(f); clear global st;

%% Another diagnostic image, looking at how well the segmentation worked...
Pthresh = 0.95;

ROIdata = roi2hist(Simg, outSeg, Pthresh);

[h, pv, ci, stats] = ttest2(ROIdata{2}, ROIdata{1});

title(sprintf('GM vs WM... T(%d) = %0.2f, p = %1.4f', stats.df, stats.tstat, pv))

set(2,'Renderer','zbuffer');
print(2,'-djpeg','-r150',...
    fullfile(localpath,['diagnostic_' aap.tasklist.main.module(aap.tasklist.currenttask.modulenumber).name '_Hist.jpg']));
try close(2); catch; end

%% Contours VIDEO of segmentations
aas_checkreg(aap,domain,indices,outSeg(1,:),Simg)
aas_checkreg(aap,domain,indices,outNSeg(1,:),tmpfile)

delete(tmpfile);
end
