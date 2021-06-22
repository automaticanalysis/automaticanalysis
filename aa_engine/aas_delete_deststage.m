function dirs2delete = aas_delete_deststage(aap,indices,dryrun)
% recursively delete destination stages

if nargin < 3, dryrun = false; end

dirs2delete = unique(aas_deststage(aap,aap.tasklist.currenttask.modulenumber,indices,{}),'stable')';
if ~dryrun
    for d = dirs2delete'
        if exist(d{1},'dir')
            aas_log(aap,false,['WARNING: deleting ' d{1}]);
            rmdir(d{1},'s'); 
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