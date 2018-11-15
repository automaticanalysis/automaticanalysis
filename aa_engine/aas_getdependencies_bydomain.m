%  Automatic analysis
%  Calculate dependencies between stages of a different type
%  Uses aap.directory_conventions.parallel_dependencies tree structure
%  Parameter 'whattoreturn':
%   'possiblestreamlocations' (default) - all places where a stream might
%   be present from the source domain that is relevant to the target domain
%   module with the specified indices
%
%   'doneflaglocations' - only specifications of modules that might produce
%   an output
%
% If 'desc' output, it is filled with a text description of the
% dependencies. This is useful for debugging
%
%  Examples:
%   pth=aas_getdependencies_bydomain(aap,'study','session',{1,2});  % source stage was at study domain, target stage is at session domain
%
% CHANGE HISTORY
%	OCT2018 [MSJ] -- added selected_sessions to hash (cf. Issue #165), cleanup
%

function [deps commonind desc]=aas_getdependencies_bydomain(aap, sourcedomain, targetdomain, indices, whattoreturn, modulenum)

if ~exist('modulenum','var')
    modulenum=[];
end
if ~exist('whattoreturn','var')
    whattoreturn='possiblestreamlocations';
end

selected_sessions = aap.acq_details.selected_sessions;
[hit resp]=aas_cache_get(aap,mfilename,sourcedomain,targetdomain,selected_sessions,indices,whattoreturn,modulenum);

if (~hit)
    % Cache miss, but actually only indices in domain tree relevant to the
    % common point of the source & target trees can affect the response
    %
    % So, find that common point and then retry the cache
    targetdomaintree=aas_dependencytree_finddomain(targetdomain,aap.directory_conventions.parallel_dependencies,{});
    sourcedomaintree=aas_dependencytree_finddomain(sourcedomain,aap.directory_conventions.parallel_dependencies,{});
    % Find the point where the source and target branches converge
    minlen=min(length(targetdomaintree),length(sourcedomaintree));
    for commonind=1:minlen
        if ~strcmp(targetdomaintree{commonind},sourcedomaintree{commonind})
            commonind=commonind-1;
            break
		end
	end
end

if (hit)
    deps=resp{1};
    commonind=resp{2};
else
    
    if length(indices)~=(length(targetdomaintree)-1)
        aas_log(aap,true,sprintf('Expected %d indicies for domain "%s" but got %d',length(targetdomaintree)-1,targetdomain,length(indices)));
	end
    
    switch(whattoreturn)
        case 'possiblestreamlocations'
            module=[];
            hasancestory=true;
            onlyreportdomain=[];
        case 'doneflaglocations'
            module=[];
            hasancestory=false;
            onlyreportdomain=sourcedomain;
        case 'doneflaglocations_thatexist'
            module=modulenum;
            hasancestory=false;
            onlyreportdomain=sourcedomain;
        otherwise
            aas_log(aap,true,sprintf('Unknown dependency calculation requested of %s',whattoreturn));
	end
    
    deps={};
    tree=aap.directory_conventions.parallel_dependencies;
    targetInd = length(targetdomaintree);
	
    % CASE 1: same branch, source module higher up
    if commonind==length(sourcedomaintree)
        if strcmp(whattoreturn,'possiblestreamlocations')
            % Bit from source to one above target
            for nodeind=commonind:targetInd
                deps{end+1}={targetdomaintree{nodeind},indices(1:nodeind-1)};
            end;
            
            % If we're looking for all possible stream locations, need to go
            % down the trees
            % Target and everything below
            for nodeind=1:targetInd
                tree=tree.(targetdomaintree{nodeind});
			end
            deps=[deps aas_dependencytree_findbranches(aap,{tree,targetdomaintree,indices},indices,onlyreportdomain,module,{},sourcedomain,hasancestory)];
        else
            deps{end+1}={sourcedomain,indices(1:length(sourcedomaintree)-1)};
		end
        % CASE 2: same branch, target module higher up
    elseif commonind==length(targetdomaintree)
        % Source and everything below
        for nodeind=1:targetInd
            tree=tree.(sourcedomaintree{nodeind});
		end
        deps=[aas_dependencytree_findbranches(aap,{tree,targetdomaintree,indices},indices,onlyreportdomain,module,{},sourcedomain,hasancestory)];
        % CASE 3: different branches
    else
        % All items below common index
        % Source and everything below
        tree=aap.directory_conventions.parallel_dependencies;
        for nodeind=1:commonind
            tree=tree.(sourcedomaintree{nodeind});
		end
        % Only pick items on or below sourcedomain
        hasancestory=false;
        deps=[aas_dependencytree_findbranches(aap,{tree,sourcedomaintree(1:commonind-1),indices(1:commonind-1)},indices(1:commonind-1),onlyreportdomain,module,{},sourcedomain,hasancestory)];
	end
    
    aas_cache_put(aap,mfilename,{deps commonind},sourcedomain,targetdomain,selected_sessions,indices,whattoreturn, modulenum);
	    
end

%% And report
if nargout==3
    desc='';
    for depind=1:length(deps)
        desc=[desc deps{depind}{1} '[ ' sprintf('%d ',deps{depind}{2}) '] '];
	end
end

end
