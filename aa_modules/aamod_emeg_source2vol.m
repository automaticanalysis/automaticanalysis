function [aap resp]=aamod_emeg_source2vol(varargin)
% Writes source contrasts to volume with 2D and 3D smoothing, assuming
% naming conventions as aamod_emeg_sourcecontrasts
%
% Danny Mitchell 02/04/08

%% check task settings, subject, block etc
[aap subblock doit resp settings]=aa_emeg_checktasksettings(mfilename('fullpath'),varargin);
if ~doit; return; end

%% find files and decide whether to run task;
files=aas_emeg_findfiles(aap,settings.InputFilter,subblock);
if isempty(files); aas_log(aap,1,sprintf('\nFound no data! (Input filter is %s)\n',settings.InputFilter)); end

for f=1:length(files);

    %% skip EEG or CDA files
    if ~isempty(regexp(files{f},'-eeg','ONCE')); continue; end
    if ~isempty(regexp(files{f},'_CDA_','ONCE')); continue; end

    fprintf('\nFile %s:',files{f});
    %% load MEG header
    try rehash; load(files{f});
    catch
        fprintf('\nFailed to load file. Will try again in 10s in case there is an access conflict...\n');
        pause(10)
        try rehash; load(files{f},'-MAT');
        catch; fprintf('\nFile may be corrupt.\n'); debugnow
        end
    end

    %% confirm (emptied?) subdirectory to save to
    [pth nam]=fileparts(D.fname);
    subdir=fullfile(D.path,nam);
    if exist(subdir,'dir')~=7; mkdir(subdir); end

    % for each inversion/coul/cosl
    %itypes={'inv','coul','cosl'};
    %itypes={'inv','coep','coip','cotp','coeJ'};
    itypes={'inv'};
    for it=1:length(itypes); % now a relic
        itype=itypes{it};
        for v=1:length(D.(itype)); % e.g. different inversions
            fprintf('\n - %s %g: %s',itype,v,D.(itype){v}.comment);
            for ctw=1:length(D.(itype){v}.contrast) % different t/f windows or evo/ind
                try
                    cname=D.(itype){v}.contrast{ctw}.cname;
                    fprintf('\n   - Contrast %g: %s',ctw,cname);
                catch aas_log(aap,true,'\nFound a source contrast without a name.\n')
                end
                contrast=D.(itype){v}.contrast{ctw};
                vname=[regexprep([D.(itype){v}.comment '_' contrast.cname],'\.','p') '.nii'];
                for e=1:length(contrast.GW) % 'contrasted' events
                    % slightly ugly output naming, but matches that from convertmat2an3D
                    edir=fullfile(subdir,sprintf('trialtype%g',e));
                    if ~exist(edir,'dir'); mkdir(edir); end
                    contrast.fname{e}=fullfile(edir, vname);
                end

                %%% minimal 2D smoothing if not MSP?! (def is 32 steps)
                if ~isempty(regexp(D.(itype){v}.comment,'MSP|GS|ARD','ONCE'))
                    contrast.SmoothingSteps2D=settings.SmoothingSteps2D;
                else contrast.SmoothingSteps2D=1;
                end
                %%%

                contrast.SmoothingMethod=settings.Smoothing3D;
                % 'fixed' to apply smoothing as normal, or 'target' to
                % estimate mean smoothness and add extra smoothing to aim
                % for target final smoothness
                % NOT YET IMPLEMENTED!!
                contrast.smooth  = settings.Smoothing3DFWHM;
                % gaussian FWHM (mm; default 8? (not 12?); Rik used 16;
                % my current default is 13; follows graph lapacian smoothing on mesh)
                contrast.display = 0;
                contrast.rms = 1;
                contrast.overwrite = settings.Overwrite;
                contrast=spm_eeg_inv_Mesh2Voxels_dm5(D,contrast);
                % add contrast name to descrip;
                % add t/f window and inverion type to file name
                % option not to overwrite
                % try to allow adjustable 3D smoothing for target smoothness - NOT YET IMPLEMENTED!!
                % option to set number of 2d smoothing steps
                % don't update D yet

                % reattach contrast (but don't save yet)
                D.(itype){v}.contrast{ctw}=contrast;

                %% mirror, then subtract or average
                if settings.mirrorsubtract
                    for e=1:length(contrast.fname)
                        Pin=contrast.fname{e};
                        [pth,nam,ext] = fileparts(Pin);
                        Lout=deblank(fullfile(pth,['l_' nam ext]));
                        Bout=deblank(fullfile(pth,['b_' nam ext]));
                        if ~exist(Lout,'file') || settings.Overwrite
                            funcstr='i1-flipdim(i1,1)';
                            fprintf('\nMirror subtraction of %s',Pin)
                            try spm_imcalc_ui(Pin,Lout,funcstr,{[],[],'float32',0});
                            catch delete(Pin);
                                error('\nFile %s seems corrupt\nIt has been deleted\n',Pin)
                            end
                        else fprintf('.')
                        end
                        if ~exist(Bout,'file') || settings.Overwrite
                            funcstr='(i1+flipdim(i1,1))/2';
                            fprintf('\nBilateral average of %s',Pin)
                            try spm_imcalc_ui(Pin,Bout,funcstr,{[],[],'float32',0});
                            catch delete(Pin);
                                error('\nFile %s seems corrupt\nIt has been deleted\n',Pin)
                            end
                        else fprintf('.')
                        end
                    end
                end

            end % next contrast type/tf window
        end % next inversion
    end % next inversion type (zul/zscoul/zucosl)

    spm_eeg_inv_save(D); % just save once to save time

end % next file

return
