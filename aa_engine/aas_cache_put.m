function [hit]=aas_cache_put(aap,name,value,varargin)
global aacache

try
    key=aas_cache_getkey(aap,varargin{:});
    aacache.(name).(key)=value;
    hit=true;
catch errcode
    hit=false;
end;