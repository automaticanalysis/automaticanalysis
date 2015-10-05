function [img, md5, inpstreamdesc] = aas_getfiles_bystream_multilevel(aap,varargin)
% Multilevel streamlocator: session/specified --> subject --> study
% Inputs:
%   aap
%   domain
%   index/indices
%   streamname
%   source: 'input'/'output' (optional, default: 'input')

source = 'input';
if strcmp(varargin{end},'input') || strcmp(varargin{end},'output')
    source = varargin{end}; 
    varargin(end) = [];
end
streamname = varargin{end};
index = varargin(1:end-1);
    
img = '';
email0 = aap.options.email; aap.options.email = 'noerror'; % silence;
[img, md5, inpstreamdesc] = aas_getfiles_bystream(aap,index{:},streamname,source);
if isempty(img)
    [img, md5, inpstreamdesc] = aas_getfiles_bystream(aap,'subject',index{2}(1),streamname,source);
end
if isempty(img)
    [img, md5, inpstreamdesc] = aas_getfiles_bystream(aap,'study',[],streamname,source);
end
aap.options.email = email0;
if isempty(img)
    aas_log(aap,1,sprintf('%s stream %s not found',source,streamname));
end
end