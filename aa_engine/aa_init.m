% Automatic analysis - initialise paths from recipe
function [aap]=aa_init(aap)

global aa

% detect running aa
if isobject(aa)
    aas_log(aap,false,'WARNING: Previous execution of aa was not closed! Killing jobs...\n')
    aas_log(aap,false,'WARNING: The path settings for both linux and MATLAB may have been modified! You may need to revise them.')
    aas_log(aap,false,'WARNING: Please, make sure that you use aa_close(aap) next time before you start a new analysis!')
    aa_close(aap,'restorewarnings','killjobs');
    aas_log(aap,false,'\nWARNING: Done!')
else
    aa = aaClass('nopath','nogreet');
end

% cleanup aaworker path
if isfield(aap.options,'aaworkercleanup') && ~isempty(aap.options.aaworkercleanup)
    aawp = aaworker_getparmpath(aap);
    for d = dir(fullfile(aawp,'aaworker*'))'
        if etime(clock,datevec(d.datenum))/(24*60*60) > aap.options.aaworkercleanup
            aas_log(aap,false,sprintf('INFO: aaworker folder %s is older than %d days...Deleting',d.name,aap.options.aaworkercleanup))
            try
                rmdir(fullfile(aawp,d.name),'s');
            catch
                aas_log(aap, false, sprintf('WARNING: Could not remove %s. Please remove manually.',fullfile(aawp,d.name)))
                pause(3)
            end
        end
    end
end

% Set UTC time function
if exist('utc_time','file')
    utc_timeFunc = @utc_time;
else
    aas_log(aap,false,'INFO: utc_time is not found. java function will be used\n')
    utc_timeFunc = @java.lang.System.currentTimeMillis;
end
aas_cache_put(aap,'utc_time',utc_timeFunc,'utils');

%% Set Paths
aas_cache_put(aap,'bcp_path',path,'system');
aas_cache_put(aap,'bcp_shellpath',getenv('PATH'),'system');

% Path for SPM
SPMDIR = '';
% - backward compatibility
if isfield(aap.directory_conventions,'spmdir') && ~isempty(aap.directory_conventions.spmdir)
    SPMDIR = aap.directory_conventions.spmdir;
end
% toolboxes
if isfield(aap.directory_conventions,'toolbox')
    tbxInd = strcmp({aap.directory_conventions.toolbox.name},'spm');
    if any(tbxInd)
        SPMDIR = aap.directory_conventions.toolbox(tbxInd).dir;
    end
end
% - path
if isempty(SPMDIR)
    if isempty(which('spm'))
        aas_log(aap,true,'You''re going to need SPM, add it to your paths manually or set in aap.directory_conventions.toolbox');
    else
        SPMDIR = spm('Dir');
    end
end
% - deployed
if isdeployed
    SPMDIR = spm('Dir');
end
% - reset
if isfield(aap.directory_conventions,'spmdir'), aap.directory_conventions.spmdir = SPMDIR; end
if isfield(aap.directory_conventions,'toolbox') && any(tbxInd)
    aap.directory_conventions.toolbox(tbxInd).name = 'spm';
    aap.directory_conventions.toolbox(tbxInd).dir = SPMDIR; 
end

% - by setting this environment variable it becomes possible to define other
%   paths relative to $SPMDIR in defaults files and task lists
setenv('SPMDIR',SPMDIR);

% - expand shell paths (before SPM so SPM can be in e.g. home directory)
aap = aas_expandpathbyvars(aap, aap.options.verbose>2);

if isfield(aap, 'spm') && isfield(aap.spm, 'defaults')
    oldspmdefaults = aap.spm.defaults;
end

SPM = spmClass(SPMDIR,'doAddToPath',true);
SPM.load;
aas_cache_put(aap,'spm',SPM);

try
    aap.spm.defaults=spm_get_defaults;
catch
    global defaults
    if isstruct(defaults)
        aas_log(aap,false,'WARNING: SPM defaults has not been found, global defaults will be used');
        aap.spm.defaults=defaults;
    else
        aap.spm.defaults = struct;
    end
end

if exist('oldspmdefaults', 'var')
    aap.spm.defaults = setstructfields(aap.spm.defaults, oldspmdefaults);
end
aap.aap_beforeuserchanges.spm.defaults = aap.spm.defaults;

% Path for matlabtools
if isfield(aap.directory_conventions,'matlabtoolsdir') && ~isempty(aap.directory_conventions.matlabtoolsdir)
    addpath(strrep(aap.directory_conventions.matlabtoolsdir,':',pathsep))
end

% Toolboxes
if isfield(aap.directory_conventions,'toolbox') && isstruct(aap.directory_conventions.toolbox)
    for TBX = reshape(aap.directory_conventions.toolbox,1,[])
        if strcmp(TBX.name,'spm'), continue; end
        aas_cache_put(aap,TBX.name,aas_inittoolbox(aap,TBX.name));
    end
end

% MNE
if isfield(aap.directory_conventions,'mnedir') && ~isempty(aap.directory_conventions.mnedir)
    if exist(fullfile(aap.directory_conventions.mnedir,'matlab'),'dir')
        addpath(aap.directory_conventions.mnedir,'matlab','toolbox');
        addpath(aap.directory_conventions.mnedir,'matlab','examples');
    end
end

% Path to GIFT
if ~isempty(aap.directory_conventions.GIFTdir)
    addpath(genpath(aap.directory_conventions.GIFTdir));
else
    % Check whether already in path, give warning if not
    if isempty(which('icatb_runAnalysis'))
       aas_log(aap,false,sprintf('GIFT not found, if you need this you should add it to the matlab path manually, or set aap.directory_conventions.GIFTdir'));
    end
end

% Path to FaceMasking
if isfield(aap.directory_conventions,'FaceMaskingdir') && ~isempty(aap.directory_conventions.FaceMaskingdir)
    addpath(genpath(fullfile(aap.directory_conventions.FaceMaskingdir,'matlab')));
else
    % Check whether already in path, give warning if not
    if isempty(which('mask_surf_auto'))
       aas_log(aap,false,sprintf('FaceMasking not found, if you need this you should add it to the matlab path manually, or set aap.directory_conventions.FaceMaskingdir'));
    end
end

% Path to LI toolbox
if isfield(aap.directory_conventions,'LIdir') && ~isempty(aap.directory_conventions.LIdir)
    addpath(aap.directory_conventions.LIdir);
else
    % Check whether already in path, give warning if not
    if isempty(which('LI'))
       aas_log(aap,false,sprintf('LI toolbox not found, if you need this you should add it to the matlab path manually, or set aap.directory_conventions.LIdir'));
    end
end


% Path to DCMTK
if isfield(aap.directory_conventions,'DCMTKdir') && ~isempty(aap.directory_conventions.DCMTKdir)
    [s,p] = aas_cache_get(aap,'bcp_shellpath','system');
    setenv('PATH',[p ':' fullfile(aap.directory_conventions.DCMTKdir,'bin')]);
end

reqpath = build_reqpath(aa.Path, aap, SPMDIR);
aas_cache_put(aap,'reqpath',reqpath,'system');

% switch off warnings
warnings(1) = warning('off','MATLAB:Completion:CorrespondingMCodeIsEmpty');
warnings(2) = warning('off','MATLAB:getframe:RequestedRectangleExceedsFigureBounds');
aas_cache_put(aap,'warnings',warnings,'system');

end

function reqpath = build_reqpath(aa_path, aap, SPMDIR)
%% Build required path list for cluster submission

% aa
reqpath = textscan(genpath(aa_path), '%s', 'delimiter', ':');
reqpath = reqpath{1};

p = textscan(path, '%s', 'delimiter', ':');
p = p{1};

% spm
reqpath = add_matched_directories(reqpath, p, SPMDIR);

% matlabtoolsdir
if isfield(aap.directory_conventions,'matlabtoolsdir') && ~isempty(aap.directory_conventions.matlabtoolsdir)
    matlabtools = textscan(aap.directory_conventions.matlabtoolsdir, '%s', 'delimiter', ':');
    matlabtools = matlabtools{1};
    for pp = matlabtools'
        if exist(pp{1},'dir')
            reqpath = [reqpath; pp{1}];
        end
    end
end

% Toolboxes
if isfield(aap.directory_conventions,'toolbox') && isstruct(aap.directory_conventions.toolbox)
    for TBX = reshape(aap.directory_conventions.toolbox,1,[])
        if isfield(TBX, 'extraparameters') ...
                && isfield(TBX.extraparameters, 'doAddToPath') ...
                && TBX.extraparameters.doAddToPath
            reqpath{end+1} = TBX.dir;
        end
    end
end

% MNE
if isfield(aap.directory_conventions,'mnedir') && ~isempty(aap.directory_conventions.mnedir)
    reqpath = add_matched_directories(reqpath, p, aap.directory_conventions.mnedir);
end

% GIFT
if ~isempty(aap.directory_conventions.GIFTdir)
    reqpath = add_matched_directories(reqpath, p, aap.directory_conventions.GIFTdir);
end


% FaceMasking
if isfield(aap.directory_conventions,'FaceMaskingdir') && ~isempty(aap.directory_conventions.FaceMaskingdir)
    reqpath = add_matched_directories(reqpath, p, aap.directory_conventions.FaceMaskingdir);
end

% LI
if isfield(aap.directory_conventions,'LIdir') && ~isempty(aap.directory_conventions.LIdir)
    reqpath = add_matched_directories(reqpath, p, aap.directory_conventions.LIdir);
end

% clean: remove empty and .git paths
reqpath = reqpath(strcmp('', reqpath)==0);
exc = cell_index(reqpath, '.git');
if exc
    reqpath(exc) = [];
end

end

function reqpath = add_matched_directories(reqpath, paths, dirname)
% Add entries from the cell array paths, whose start matches dirname, to the end of the reqpath cell array.

% In paths, the entries do not have any trailing filesep
% Make sure dirname itself is also a match
if strcmp(dirname(end), filesep)
    % Remove trailing filesep
    dirname = dirname(1:end-1);
end

p_ind = cell_index(paths, dirname);
nr_add = length(p_ind);
reqpath(end+1:end+nr_add) = paths(p_ind);

end
