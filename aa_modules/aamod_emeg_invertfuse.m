function [aap resp]=aamod_emeg_invertfuse(varargin)
% Fuse inversions (or group inversions) that have been estimated seperately
% for different sensor types. Parallelised over sessions of each subject.
%
% Danny Mitchell 02/04/08; 07/11/08
% Adapted from Rik's cbu_meeg_spm5_pipeline.m
%
% To do:
% Needs to be generalised to accept fusion of EEG
% Gereralise to fuse different inversions in same file
% Currently uses same inversion setting as were specified for group
% inversion of seperate sensor type; potentially may want to be changed...
% Allow overwriting/non-overwriting?

%% check task settings, subject, block etc
[aap subblock doit resp settings]=aa_emeg_checktasksettings(mfilename('fullpath'),varargin);
if ~doit; return; end

typ={'-mags','-grds'};
Ntyp=length(typ);

% collect file names
for typn=1:Ntyp
    filt=sprintf('^sace.*_BLOCK.*%s_grp\\.mat$',typ{typn});
    files=aas_emeg_findfiles(aap,filt,subblock);
    for fstem=1:length(files)
        fnam{fstem,typn}=files{fstem};
    end
end
missing=cellfun(@isempty,fnam(:));
if any(missing)
    fprintf('\n *** Failed to fuse inversions from multiple sensor types. ***\nFound following files:\n')
    disp(char(fnam(~missing)));
    aas_log(aap,1,'Fusion inversion failed because could not match up all sensor types!');
end

for fstem=1:size(fnam,1)
    
    clear DD;
    ids={};
    for typn=1:Ntyp
        DD{typn} = spm_eeg_ldata(fnam{fstem,typn});
        for v=1:length(DD{typn}.inv)
            ids=[ids; DD{typn}.inv{v}.comment];
            DD{typn}.inv{v}.inverse.Nm = [];   % Reset any parameters from last inversion
            DD{typn}.inv{v}.inverse.Nr = [];   % Reset any parameters from last inversion
            % DD{typn}.inv{v}.inverse.QP =[]; (Important not to reset the source priors from the grp
            % inversion above)
        end
    end
    ids=unique(ids);
    
    for idind=1:length(ids) % fuse all inversion types found in file
        for typn=1:Ntyp
            % find index of this inversion type for each sensor file
            status='unknown';
            for v=1:length(DD{typn}.inv)
                try
                    if strcmp(DD{typn}.inv{v}.comment,ids{idind});
                        status='found';break
                    end
                catch
                end
            end
            if strcmp(status,'unknown');
                aas_log(aap,1,'\n Failed to find inversion with ID: %s\n in file: %s',id,fnam{fstem,typn});
            else
                DD{typn}.val=v;
            end            
        end

    D = spm_eeg_invert_fuse(DD);  % note necessary so not overwritten by next val
    % Inversion appended and saved; adds suffix '_fused' to first filename
    
    end % next inversion type

end % next file 

% %%%%%%%%%%%%%%%%% Quick look at mean over subjects...
%
% clear S; S.con = [1 -1]; S.boi{1} = [140 190]; S.boi{2} = []; % boi{2}=[] is all freqs
% S.fig = 6;
% for val=1:Nval
%     S.val = val;
%     for typn=1:(Ntyp+fusionflag)
%         S.fig = S.fig+1;
%         S.P = strvcat(grpfinalnam{typn,:});
%         meg_inv_display_con(S);
%         drawnow
%     end
% end

