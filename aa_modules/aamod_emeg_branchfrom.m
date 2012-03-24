function [aap resp]=aamod_emeg_branchfrom(varargin)
% Copy files from old analysis
% Danny Mitchell 13/11/08

%% check task settings, subject, block etc
[aap subblock doit resp settings]=aa_emeg_checktasksettings(mfilename('fullpath'),varargin);
if ~doit; return; end

%% group level
fprintf('\nAnalysis %s > %s', settings.OldAnalysis, aap.acq_details.root)

if exist(settings.OldAnalysis,'dir')~=7
    pth=fileparts(regexprep(aap.acq_details.root,'/$',''));
    settings.OldAnalysis=fullfile(pth,settings.OldAnalysis);
    if exist(settings.OldAnalysis,'dir')~=7
        aas_log(aap,1,sprintf('\nFailed to find old analysis: %s\n', settings.OldAnalysis))
    end
end

if settings.CopyEvents
    fprintf('\n - Copying events folder')
    if settings.Overwrite
        copyfile(fullfile(settings.OldAnalysis,'events'),fullfile(aap.acq_details.root,'events'),'f');
    else copyfile(fullfile(settings.OldAnalysis,'events'),fullfile(aap.acq_details.root,'events'));
    end
end

if settings.CopyBadChans && exist(fullfile(settings.OldAnalysis,'bad_chans.xls'),'file')
    if settings.Overwrite
        copyfile(fullfile(settings.OldAnalysis,'bad_chans.xls'),aap.acq_details.root,'f');
    else copyfile(fullfile(settings.OldAnalysis,'bad_chans.xls'),aap.acq_details.root);
    end
end

%% subject level
for s=1:length(aap.acq_details.subjects)
    subname=aap.acq_details.subjects(s).megname;
    newsubdir=fullfile(aap.acq_details.root,subname);
    fprintf('\n - Subject %s',subname)

    if settings.CopyEvents
        fprintf('\n   - Copying events folder')
        if settings.Overwrite
            copyfile(fullfile(settings.OldAnalysis,subname,'events'),fullfile(newsubdir,'events'),'f');
        else copyfile(fullfile(settings.OldAnalysis,subname,'events'),fullfile(newsubdir,'events'));
        end
    end

    if settings.CopyStructurals
        fprintf('\n   - Copying structurals folder')
        if settings.Overwrite || ~exist(fullfile(newsubdir,'structurals'),'dir')
            copyfile(fullfile(settings.OldAnalysis,subname,'structurals'),fullfile(newsubdir,'structurals'),'f');
        end
    end

%% session level
    for b=aap.acq_details.selected_sessions
        filt=strrep(settings.InputFilter,'BLOCK','');
        block=aap.acq_details.sessions(b).name;
        fprintf('\n   - Block %s',block)
        newblockdir=fullfile(newsubdir,block);
        oldblockdir=fullfile(settings.OldAnalysis,subname,block);
        if settings.CopyEvents
            if settings.Overwrite
                copyfile(fullfile(oldblockdir,'events'),fullfile(newblockdir,'events'),'f');
            else copyfile(fullfile(oldblockdir,'events'),fullfile(newblockdir,'events'));
            end
        end

        files=cellstr(spm_select('List',oldblockdir,filt));
        for f=1:length(files)
            original=fullfile(oldblockdir,files{f});
            copied=fullfile(newblockdir,files{f});
            if ~exist(copied,'file') || settings.Overwrite==1
                copyfile(original,copied,'f')
                
                % check and update D.path if necessary
                try
                    rehash
                    clear D
                    load(copied)
                    if isfield(D,'path')
                        D.path=fileparts(copied);
                        save(copied,'D');
                    end
                catch
                end
                
            end
            fprintf('.')
        end

    end % next block
end % next subject

return