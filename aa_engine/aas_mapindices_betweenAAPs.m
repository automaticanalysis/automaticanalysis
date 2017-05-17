function map = aas_mapindices_betweenAAPs(localAA, remoteAA)
% function map = aas_mapindices_betweenAAs(curDomain, localAA, remoteAA)
%
% This function creates a mapping between the indices of different AAP
% structures. For example, the local AAP might contain a subset of subjects
% and sessions of the remote AAP, and a result, subject index '1' will access
% different data depending on which AAP you are referring to. 
%
% The input 'curDomain' determines the point in the dependency tree from
% whcih we descend. If you provide 'study', you will get mappings for
% everything. If you specify 'study' you will get the mappings for all
% domain contained in a subject, as well as the mapping for the subject
% domain.
%
% The returned map is a structure with a field for every domain from the
% requested domain on down. Each field is a vector, that when indexed with
% local index, returns the index in the remote AA.
%
% For example, if the remote AA has three subjects: 'SubA' 'SubB' 'SubC',
% and the local AA has only two of those: 'SubC' 'SubA',
% but all the subjects have three sessions in the correct order,
% then the map looks like this:
%
%       study: 1
%     subject: [3 1]
%     session: [1 2 3]
%        ...
% _______________________________________________________________________
%
% created by cwild 2014-09-09
%
%

% I guess this isn't really needed, but might be useful! I don't remember
% why I started it this way...
[domainFound, localTree] = recursiveDescent(localAA.directory_conventions.parallel_dependencies, 'study');

if ~domainFound
    aas_log(localAA, 1, sprintf('Invalid requested domain ''%s'': not present in aap.directory_conventions.parallel_dependencies', requestedDomain));
end

map = recursiveMap('study', localTree, localAA, remoteAA, struct());

end

function [foundIt, subTree] = recursiveDescent(localTree, domainName)

if isstruct(localTree)
    subTrees = fieldnames(localTree);
    if ismember(domainName, subTrees)
        subTree = localTree.(domainName);
        foundIt = true;
    else
        for i = 1 : numel(subTrees)
            [foundIt, subTree] = recursiveDescent(localTree.(subTrees{i}), domainName);
            if(foundIt), return;
            end
        end
    end
else
    foundIt = false;
    subTree = [];
end

end

function map = recursiveMap(curDomain, localTree, localAA, remoteAA, map)

remoteDep = aas_dependencytree_finddomain(curDomain, remoteAA.directory_conventions.parallel_dependencies, {});
if isempty(remoteDep), return; end

localDomainNames = aas_getNames_bydomain(localAA, curDomain);
remoteDomainNames = aas_getNames_bydomain(remoteAA, curDomain);

if ~strcmp(curDomain, 'study')
    [junk,  map.(curDomain)] = ismember(localDomainNames{1}, remoteDomainNames{1});
    
else
    map.(curDomain) = 1;
end

if isstruct(localTree)
    subTrees = fieldnames(localTree);
    for tI = 1 : numel(subTrees)
        map = recursiveMap(subTrees{tI}, localTree.(subTrees{tI}), localAA, remoteAA, map);
    end
end

end

