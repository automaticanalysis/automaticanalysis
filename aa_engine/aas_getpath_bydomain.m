% Automatic analysis - get path using domain (subject, session etc)
%  Uses aap.directory_conventions.parallel_dependencies tree structure
%
% examples:
%  pth=aas_getpath_bydomain(aap,'session',[1,2]);  % subject 1, session 2
%  pth=aas_getpath_bydomain(aap,'session',[1,2],4);  % subject 1, session 2 for module 4
%  pth=aas_getpath_bydomain(aap,'searchlight',[1,2,3],'s3');   % path of searchlight on s3
%  pth=aas_getpath_bydomain(aap,'session',[1,2],'s3',4);   % path on s3 for module 4

function [pth]=aas_getpath_bydomain(aap,domain,indices,varargin)


domaintree=aas_dependencytree_finddomain(domain,aap.directory_conventions.parallel_dependencies,{});

if length(indices)~=(length(domaintree)-1)
    aas_log(aap,true,sprintf('Expected %d indicies for domain "%s" but got %d',length(domaintree)-1,domain,length(indices)));
end;

pth=aas_getstudypath(aap,varargin{:});

for ind=1:(length(domaintree)-1)
    pth=fullfile(pth,aas_getdirectory_bydomain(aap,domaintree{ind+1},indices(ind)));
end;



end
