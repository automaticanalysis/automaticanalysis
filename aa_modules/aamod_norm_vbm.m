% AA module - normalisation using normalise or two pass procedure with segment
% [aap,resp]=aamod_norm_vbm(aap,task,i)
% Depending on aap.tasksettings.aamod_norm_noss.usesegmentnotnormalise
% i=subject num
% Rhodri Cusack & Daniel Mitchell MRC CBU 2006
% based on originals by Rik Henson, Matthew Brett
% @@@ THIS IS NOT YET TRANSFORMED TO AA4 @@@

function [aap,resp]=aamod_norm_vbm(aap,task,i)

resp='';

switch task
    case 'report'

    case 'doit'

        global defaults;
        defs =defaults.normalise;
        defs.estimate.weight = '';

        % Template image; here template image not skull stripped
        temp_imgs = aap.directory_conventions.T1template;

        clear imgs;
        subj_dir = aas_getsubjpath(aap,i);
        % Get structural for this subject
        structdir=fullfile(subj_dir,aap.directory_conventions.structdirname);
        %structfilter=fullfile(structdir,['s aap.acq_details.subjects(i).structuralfn '*.nii']); % doesn't work for sMR* images
        structfilter=fullfile(structdir,'s*.nii');
        subj.P = dir(structfilter);
        
        if (length(subj.P)>1)
            aas_log(aap,0,sprintf('Found more than one structural so using first:\n%s',subj.P));
        end
        
        if isempty(subj.P) || strcmp(subj.P, '/')
            aas_log(aap,1,sprintf('No structural found with filter %s',structfilter));
        end

        % Set the mask for subject to empty by default
        subj.objmask = ''; % object mask

        % Get the images we are going to reslice
        subj.PP = subj.P;

        % log which SPM function we are using
        aas_log(aap,0,sprintf('Running aamod_norm_vbm with %s', which('spm_preproc')));
        
        
        % call the SPM segment or normalize function to do the work
        % 2 stage process, as proposed by RH, to increased robustness [djm 13/03/06]
        %%%%%%%% 1st pass:
        estopts.regtype='';    % turn off affine:
        out = spm_preproc(fullfile(structdir,subj.P(1).name),estopts);
        [sn,isn]   = spm_prep2sn(out);

        % only write out attenuation corrected image
        writeopts.biascor = 1;
        writeopts.GM  = [0 0 0];
        writeopts.WM  = [0 0 0];
        writeopts.CSF = [0 0 0];
        writeopts.cleanup = [0];
        spm_preproc_write(sn,writeopts);

        %%%%%%%% 2nd pass using attenuation corrected image
        % mstruc should be the attenuation corrected image

        % look for m prefixed filename
        mstruc.P = dir( fullfile(structdir,'ms*.nii'));
        if (length(mstruc.P)>1)
            aas_log(aap,0,sprintf('Found more than one attenuated structural so using first:\n%s',subj.P));
        end;

        estopts.regtype='mni';    % turn on affine again
        estopts.samp = aap.tasksettings.aamod_norm_vbm.samp;
        out = spm_preproc(fullfile(structdir,mstruc.P(1).name),estopts);
        [sn,isn]   = spm_prep2sn(out);

        % write out GM and WM native + unmod normalised
        writeopts.biascor = 1;
        writeopts.GM  = [1 1 1];	% assume GM(2) means unmod
        writeopts.WM  = [1 1 1];
        writeopts.CSF = [1 0 0];
        spm_preproc_write(sn,writeopts);

        subj.matname = fullfile(structdir,[spm_str_manip(mstruc.P(1).name,'sd') '_seg_sn.mat']);
        subj.invmatname = fullfile(structdir,[spm_str_manip(mstruc.P(1).name,'sd') '_seg_inv_sn.mat']);
        savefields(subj.matname,sn);
        savefields(subj.invmatname,isn);

        spm_write_sn(fullfile(structdir,subj.PP(1).name),subj.matname,defs.write);


        % Now save graphical check
        figure(spm_figure('FindWin', 'Graphics'));
        % added graphical check for when segment is used [djm 20/01/06]
        myvols=spm_vol(char(aap.directory_conventions.T1template, ... % template T1
            fullfile(structdir,subj.P(1).name), ... % native T1
            fullfile(structdir,strcat('w',subj.P(1).name)), ... % normalised T1
            fullfile(structdir,strcat('c1m',subj.P(1).name)))); % native grey matter segmentation
        spm_check_registration(myvols)

        ann1=annotation('textbox',[.05 .96 .9 .03],'HorizontalAlignment','center','Color','r','String',strcat('Subject:...',subj.P(1).name,',  processed on:...',date));
        ann2=annotation('textbox',[.1 .891 .3 .025],'HorizontalAlignment','center','Color','r','String','T1 template');
        ann3=annotation('textbox',[.6 .891 .3 .025],'HorizontalAlignment','center','Color','r','String','Native T1');
        ann4=annotation('textbox',[.1 .413 .3 .025],'HorizontalAlignment','center','Color','r','String','Normalised T1');
        ann5=annotation('textbox',[.6 .413 .3 .025],'HorizontalAlignment','center','Color','r','String','Native segmented grey matter');
        print('-djpeg',fullfile(subj_dir,'diagnostic_aamod_norm_noss'));

        print('-djpeg',fullfile(subj_dir,'diagnostic_aamod_norm_noss'));
       
        delete(ann1); delete(ann2); delete(ann3);delete(ann4);delete(ann5);
       
end;

%------------------------------------------------------------------------
function savefields(fnam,p)
if length(p)>1, error('Can''t save fields.'); end;
fn = fieldnames(p);
if numel(fn)==0, return; end;
for i=1:length(fn),
    eval([fn{i} '= p.' fn{i} ';']);
end;
if str2double(version('-release'))>=14,
    save(fnam,'-V6',fn{:});
else
    save(fnam,fn{:});
end;

return;
%------------------------------------------------------------------------










