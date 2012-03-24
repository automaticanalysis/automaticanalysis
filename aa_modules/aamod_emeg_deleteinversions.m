function [aap resp]=aamod_emeg_deleteinversions(varargin)

% Danny Mitchell 02/04/08; 07/11/08

%% check task settings, subject, block etc
[aap subblock doit resp settings]=aa_emeg_checktasksettings(mfilename('fullpath'),varargin);
if ~doit; return; end
try settings=rmfield(settings,'specialrequirements'); catch; end

%% find files and decide whether to run task;
files=aas_emeg_findfiles(aap,settings.InputFilter,subblock);
if isempty(files); aas_log(aap,1,sprintf('\nFound no data! (Input filter is %s)\n',settings.InputFilter)); end
clear varargin

%% for each file 
for f=1:length(files);

    load(files{f});
    try
        fprintf('\n%s: deleting %g inversions',files{f},length(D.inv))
        if ~isempty(D.inv)
            D.inv=struct([]);
            if isfield(D,'coul'); D=rmfield(D,'coul'); end
            if isfield(D,'cosl'); D=rmfield(D,'cosl'); end
            save(fullfile(D.path,D.fname),'D');
        end
    catch
        fprintf('\n%s: no inversions',files{f})
    end
    
end % next file

