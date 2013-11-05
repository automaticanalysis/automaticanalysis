% Automatic analysis
%  Caculatuate dependencies between stages of a different type
%  Uses aap.directory_conventions.parallel_dependencies tree structure
%  Parameter 'whattoreturn':
%   'possiblestreamlocations' (default) - all places where a stream might
%   be present from the source domain that is relevant to the target domain
%   module with the specified indices
%
%   'doneflaglocations' - only specifications of modules that might produce
%   an output
%
%  Examples:
%   pth=aas_getdependencies_bydomain(aap,'study','session',{1,2});  % source stage was at study domain, target stage is at session domain
%
%
function [deps commonind]=aas_getdependencies_bydomain(aap, sourcedomain, targetdomain, indices, whattoreturn, modulenum)

if ~exist('modulenum','var')
    modulenum=[];
end;
if ~exist('whattoreturn','var')
    whattoreturn='possiblestreamlocations';
end;


[hit resp]=aas_cache_get(aap,mfilename,sourcedomain,targetdomain,indices,whattoreturn,modulenum);
hit=false;

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
        end;
    end;
    
    
    
    if (commonind<length(targetdomaintree))
        cacheindices=indices(1:(commonind)); % CW hacky fix - indices(1:(commonind-1));
        [hit resp]=aas_cache_get(aap,mfilename,sourcedomain,targetdomaintree{commonind},cacheindices,whattoreturn,modulenum);
    else
        cacheindices=indices;
    end;
end;

if (hit)
    deps=resp{1};
    commonind=resp{2};
else
    
    if length(indices)~=(length(targetdomaintree)-1)
        aas_log(aap,true,sprintf('Expected %d indicies for domain "%s" but got %d',length(sourcedomaintree)-1,targetdomain,length(indices)));
    end;
    
    
    switch(whattoreturn)
        case 'possiblestreamlocations'
            targetInd = length(targetdomaintree);
            % CASE 1: same branch, source module higher up
            if commonind==length(sourcedomaintree)
                % Bit from source to one above target
                for nodeind=commonind:targetInd-1
                    deps{end+1}=targetdomaintree(1:nodind);
                end;
                % Target and everything below
                for nodeind=1:targetInd
                    tree=tree.(targetdomaintree{nodeind});
                end;  
                deps=[deps aas_dependencytree_findbranches(aap,{tree,targetdomaintree(2:end),indices(targetInd:end)},indices(1:(targetInd-1)))]; 
            % CASE 2: same branch, target module higher up
            elseif commonind==length(targetdomaintree)
                % Source and everything below
                tree=aap.directory_conventions.parallel_dependencies;
                for nodeind=1:targetInd
                    tree=tree.(sourcedomaintree{nodeind});
                end;  
                deps=[aas_dependencytree_findbranches(aap,{tree,sourcedomaintree(2:end),indices(targetInd:end)},indices(1:(targetInd-1)))]; 
            % CASE 3: different branches 
            else
                % All items below common index
                % Source and everything below
                tree=aap.directory_conventions.parallel_dependencies;
                for nodeind=1:commonind
                    tree=tree.(targetdomaintree{nodeind});
                end;
                deps=[aas_dependencytree_findbranches(aap,{tree,sourcedomaintree(2:commonind),indices(1:(commonind-1))},indices(1:(targetInd-1)))];
            end;
%            indices(commonind:end)},indices(1:(commonind-1))
        case 'doneflaglocations_thatexist'
            % Is the common point in the source & target trees the one we're
            % looking for?
            if strcmp(sourcedomain,sourcedomaintree{commonind})
                deps={{sourcedomain indices(1:(commonind-1))}};
            else
                % No, so lets go down the tree and look for it, providing every
                % possible sub index combination
                tree=aap.directory_conventions.parallel_dependencies;
                for targetind=1:commonind
                    tree=tree.(sourcedomaintree{targetind});
                end;
                deps=[aas_dependencytree_findbranches(aap,tree,indices(1:(commonind-1)),sourcedomain,modulenum)];
            end;
        case 'doneflaglocations'
            
            % Is the common point in the source & target trees the one we're
            % looking for?
            if strcmp(sourcedomain,sourcedomaintree{commonind})
                deps={{sourcedomain indices(1:(commonind-1))}};
            else
                % No, so lets go down the tree and look for it, providing every
                % possible sub index combination
                tree=aap.directory_conventions.parallel_dependencies;
                for targetind=1:commonind
                    tree=tree.(sourcedomaintree{targetind});
                end;
                deps=[ aas_dependencytree_findbranches(aap,tree,indices(1:(commonind-1)),sourcedomain)];
            end;
        otherwise
            aas_log(aap,true,sprintf('Unknown dependency calculation requested of %s',whattoreturn));
    end;
    
    aas_cache_put(aap,mfilename,{deps commonind},sourcedomain,targetdomain,indices,whattoreturn,modulenum);
    % also store in the cache indexed by the truncated portion
    aas_cache_put(aap,mfilename,{deps commonind},sourcedomain,targetdomaintree{commonind},cacheindices,whattoreturn,modulenum);
    
end







