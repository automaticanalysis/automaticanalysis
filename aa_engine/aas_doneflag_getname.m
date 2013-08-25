function [doneflagname]=aas_doneflag_getname(aap,stage)
% allow full path of module to be provided

doneflagname=['done_' aas_getstagetag(aap,stage)];
