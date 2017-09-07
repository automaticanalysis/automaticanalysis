% function [hit resp]=aas_cache_get(aap,name,varargin)
% Get item from cache using arguments in name and varargin to generate key
% name must be suitable as a matlab structure field name 

function [hit resp]=aas_cache_get(aap,name,varargin)
global aacache

try
    key=aas_cache_getkey(aap,varargin{:});
    resp=aacache.(name).(key);
    hit=true;
catch errcode
    hit=false;
    resp='';
end;