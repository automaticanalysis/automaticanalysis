function aa_close
% cleanup 

global aacache;

if isstruct(aacache) && isfield(aacache,'bcp_path')
    path(aacache.bcp_path);
end

clear global aacache aaparallel aaworker taskqueue localtaskqueue;

end

