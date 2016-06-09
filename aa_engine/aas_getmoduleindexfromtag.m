function [index, name, id] = aas_getmoduleindexfromtag(aap, stagetag)
%  From a stage tag (e.g., 'aamod_realign_00001'), returns the index of
%  that module in aap.tasklist.main.module, as well as the non-qualified
%  name and the instance of that type of module (i.e., the _00001 part).
%
%  Stages are processing stages (e.g., aamod_realign)
%  From aa_ver_3 onwards, it is possible to have multiple instances of a
%  single module in a task list - so realign could occur twice in theory
%  The done flags and dependencies for the mulitple instances are
%  differentiated by the addtion of _0002 _0003 etc
%  The stage tag is the stage name, with this index (e.g.,
%  aamod_realign_0002)
%
%  RC 15/2/2010
%  From version 4.0, even when only one stage, index is added, as this way
%  if a later stage is added that repeats one already present, previous one
%  will be identified

m = regexp(stagetag, '_\d{5,5}');
if isempty(m), aas_log(aap,true,sprintf('Invalid stage tag: %s', stagetag)); end
    
name = stagetag(1:m-1);
[path name] = fileparts(name); % allow fully qualified name, trim the path 

try
    id = str2num(stagetag(m+1:end));
catch
    aas_log(aap,true,sprintf('Invalid stage tag: %s', stagetag));
end

nameMatch = find(strcmp(name, {aap.tasklist.main.module.name}));
if isempty(nameMatch)
    aas_log(aap,true,sprintf('%s not found in aap.tasklist.main', name));
end

idMatch = find(id == [aap.tasklist.main.module(nameMatch).index] );
if isempty(idMatch)
    aas_log(aap,true,sprintf('%s with index=%d not found aap.tasklist.main', name, id));
elseif (length(idMatch)>1)
    aas_log(aap,true,sprintf('Somehow we have multiple %s with index=%d. That shouldn''t happen, so something is probably wrong with aarecipe.', name, id));
end

index = nameMatch(idMatch);
if index > length(aap.tasklist.main.module)
    aas_log(aap,true,sprintf('Somehow we have a module index that is greater than the number of modules in aap.tasklist.main'));
end


end


