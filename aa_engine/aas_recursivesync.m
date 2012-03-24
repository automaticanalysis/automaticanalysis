%% function [aap todo]=aas_recursivesync(aap,srcroot,destroot,deletefromdest)
%   deletefromdest - if true, delete files in dest tree not present in src

function [aap todo]=aas_recursivesync(aap,srcroot,destroot,deletefromdest)

[aap src]=aas_recursivedir(aap,srcroot);
[aap dest]=aas_recursivedir(aap,destroot);

if (~exist('deletefromdest','var'))
    deletefromdest=false;
end;

todo.todeletefiles={};
todo.todeletedirs={};
todo.filestocopy={};

if isstruct(dest)
    todo.dirstomake={};
else
    todo.dirstomake={''};
end;

[aap todo]=aas_recursivesync_fromstruct(aap,todo,srcroot,destroot,'',src,dest);

% Now do sync
%  delete inappropriate files in dest
if (deletefromdest)
    % delete files
    for i=1:length(todo.todeletefiles)
        cmd=['rm ' fullfile(destroot,todo.todeletefiles{i})];
        aas_shell(cmd);
    end
    
    % delete dirs, no such thing in s3
    for i=1:length(todo.todeletedirs)
        cmd=['rmdir ' fullfile(destroot,todo.todeletedirs{i})];
        aas_shell(cmd);
    end
end;

%  and make directories (no such thing in s3)
for i=1:length(todo.dirstomake)
    aas_makedir(aap,fullfile(destroot,todo.dirstomake{i}));
end

%  and copy (for s3, use bulk upload command)
for i=1:length(todo.filestocopy)
    fle=todo.filestocopy{i};
    cmd=['cp ' fullfile(srcroot,fle) ' ' fullfile(destroot,fle)];
    aas_shell(cmd);
end

end


%% And recursively calculate what needs to be done

function [aap todo]=aas_recursivesync_fromstruct(aap,todo,srcroot,destroot,relpath,src,dest)

for ind=1:length(src)
    % find corresponding item in dest
    if (~isempty(dest))
        destpos=strmatch(src(ind).name,{dest.name});
    else
        destpos=[];
    end;
    
    % mark this item as having a match and hence not for deletion
    if (~isempty(destpos) && dest(destpos).isdir==src(ind).isdir)
        dest(destpos).gotmatch=true;
    else
        destpos=[];
    end;
    
    if (src(ind).isdir)
        if (isempty(destpos))
            todo.dirstomake{end+1}=fullfile(relpath,src(ind).name);
        end;
        if (~isempty(destpos))
            [aap todo]=aas_recursivesync_fromstruct(aap,todo,srcroot,destroot,fullfile(relpath,src(ind).name),src(ind).subdir,dest(destpos).subdir);
        else
            [aap todo]=aas_recursivesync_fromstruct(aap,todo,srcroot,destroot,fullfile(relpath,src(ind).name),src(ind).subdir,[]);
        end;
    else
        if (isempty(dest) || isempty(destpos) || src(ind).datenum>dest(destpos).datenum || src(ind).bytes~=dest(destpos).bytes)
            todo.filestocopy{end+1}=fullfile(relpath,src(ind).name);
        end;
    end;
end;

% Now delete everything in dest that doesn't match src
for ind=1:length(dest)
    if (~isfield(dest(ind),'gotmatch') || isempty(dest(ind).gotmatch))
        if (dest(ind).isdir)
            todo.todeletedirs{end+1}=fullfile(relpath,dest(ind).name);
        else
            todo.todeletefiles{end+1}=fullfile(relpath,dest(ind).name);
        end;
    end;
end;
end

