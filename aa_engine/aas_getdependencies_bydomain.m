% Automatic analysis - get path using domain (subject, session etc)
%  Uses aap.directory_conventions.parallel_dependencies tree structure
%
% examples:
%  pth=aas_getdependencies_bydomain(aap,'study','session',{1,2});  % source stage was at study domain, target stage is at session domain

function [deps]=aas_getdependencies_bydomain(aap, sourcedomain, targetdomain, indices, varargin)

targetdomaintree=aas_dependencytree_finddomain(targetdomain,aap.directory_conventions.parallel_dependencies,{});
sourcedomaintree=aas_dependencytree_finddomain(sourcedomain,aap.directory_conventions.parallel_dependencies,{});

if length(indices)~=(length(targetdomaintree)-1)
    aas_log(aap,true,sprintf('Expected %d indicies for domain "%s" but got %d',length(sourcedomaintree)-1,domain,length(indices)));
end;

% Where the trees are the same at the start, we're only dependent on the
% last common item

minlen=min(length(targetdomaintree),length(sourcedomaintree));
for commonind=1:minlen
    if ~strcmp(targetdomaintree{commonind},sourcedomaintree{commonind})
        break
    end;
end;

% First, lets go up the tree as far as the source domain
deps={};
for targetind=length(targetdomaintree):-1:commonind
    deps{end+1}={targetdomaintree{targetind}, indices(1:(targetind-1))};
end;

% And down the tree, sucking everything along the way
tree=aap.directory_conventions.parallel_dependencies;
for targetind=1:length(targetdomaintree)
    tree=tree.(targetdomaintree{targetind});
end;

% 
deps=[deps aas_dependencytree_findbranches(aap,tree,deps,indices)];

end







