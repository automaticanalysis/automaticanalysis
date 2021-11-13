function dirs2delete = aa_delete_deststage(aap,indices,varargin)
% recursively delete/backup destination stages
%
% It removes or backs up existing (intermediate) steps of a workflow before
% the re-execution of these steps (with different settings).
%
% Compulosory positional arguments:
%   aap:        aap structure actualised for the stage whose destinations 
%               are to be removed/backed-up
%   indices:    array of indices corresponding to the domain of the stage,
%               e.g.
%                   study   - []
%                   subject - [subj]
%                   session - [subj sess]
% 
% Optional parameters (key-value pairs):
%   includeCurrent: whether to remove/back up the current stage, as well
%                   (default: false)
%   dryrun:         whether to only compile the list of stage folders to be
%                   removed/backed-up (default: false)
%   backupFolder:   path to the folder in which the stages are backed-up,
%                   The folder will be created if does not exist. If empty,
%                   the operation is delete rather than back-up. (default: '')

argParse = inputParser;
argParse.addParameter('includeCurrent',false,@(x) islogical(x) || isnumeric(x));
argParse.addParameter('dryrun',false,@(x) islogical(x) || isnumeric(x));
argParse.addParameter('backupFolder','',@ischar);
argParse.parse(varargin{:});

if ~isempty(argParse.Results.backupFolder), aas_makedir(aap,argParse.Results.backupFolder); end

dirs2delete = unique(aas_deststage(aap,aap.tasklist.currenttask.modulenumber,indices,{}),'stable')';
if argParse.Results.includeCurrent
    dirs2delete = [cellstr(aas_getpath_bydomain(aap,aap.tasklist.currenttask.domain,indices,aap.tasklist.currenttask.modulenumber));dirs2delete];
end

if ~argParse.Results.dryrun
    for d = dirs2delete'
        if exist(d{1},'dir')
            if ~isempty(argParse.Results.backupFolder)
                dest = strrep(d{1},aap.internal.aap_initial.acq_details.root,argParse.Results.backupFolder);
                aas_log(aap,false,['INFO: moving ' d{1} ' to ' argParse.Results.backupFolder]);
                movefile(d{1},dest);
            else
                aas_log(aap,false,['WARNING: deleting ' d{1}]);
                rmdir(d{1},'s');
            end
        end
    end
end
end

function dirlist = aas_deststage(aap,k,indices,dirlist)
srcdomain = aap.tasklist.currenttask.domain;
for o=1:numel(aap.internal.outputstreamdestinations{k}.stream)
    dest_aap = aas_setcurrenttask(aap,aap.internal.outputstreamdestinations{k}.stream(o).destnumber);
    destdomain = dest_aap.tasklist.currenttask.domain;
    if is_ancestordomain(aap,destdomain,srcdomain)
        dep = aas_getdependencies_bydomain(aap,destdomain,srcdomain,indices);
        dirlist{end+1} = aas_getpath_bydomain(dest_aap,dep{1}{1},dep{1}{2});
    else
        dirlist{end+1} = aas_getpath_bydomain(dest_aap,srcdomain,indices);
    end
    
    dirlist = aas_deststage(aap,...
        aap.internal.outputstreamdestinations{k}.stream(o).destnumber,...
        indices,dirlist);
end
end

function val = is_ancestordomain(aap,domain1,domain2)
tree = get_depenencytreefrom(domain1,aap.directory_conventions.parallel_dependencies);
val = ~isempty(get_depenencytreefrom(domain2,tree));
end

function tree = get_depenencytreefrom(domain,trunk)
tree = [];
if ~isstruct(trunk), return; end

fields = fieldnames(trunk);
if any(strcmp(fields,domain)), tree = trunk.(domain);
else
    for f = fields'
        tree = get_depenencytreefrom(domain,trunk.(f{1}));
        if ~isempty(tree), break; end
    end
end
end