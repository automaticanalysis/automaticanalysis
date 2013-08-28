function [md5_hex]=aas_cache_getkey(aap,varargin)

key=[];
md=java.security.MessageDigest.getInstance('MD5');
for arg=1:length(varargin)
    
    if isnumeric(varargin{arg})
        field=['n' sprintf('%f',varargin{arg})];
        field(field==' ')='_';
    elseif ischar(varargin{arg})
        field=['c' varargin{arg}];
    else
        aas_log(aap,false,sprintf('Caching does not support data type %s',class(varargin(arg))));
        return;
    end;
    md.update(uint8(field),0,length(field));
end;
md5_hex=sprintf('%02x',uint8(md.digest));
md5_hex=['k' md5_hex];