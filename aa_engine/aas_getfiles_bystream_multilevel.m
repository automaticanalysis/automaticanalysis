function [img, md5, streamdesc, domain] = aas_getfiles_bystream_multilevel(aap,varargin)
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

domains = {index{1}, aas_getsesstype(aap), 'subject', 'study'};
[junk, ind] = unique(domains,'legacy');
domains = domains(sort(ind));
ind0 = cell_index(domains,index{1}); ind0 = min(ind0);
if ind0, domains = domains(ind0:end); end

img = ''; d = 0;
verbose0 = aap.options.verbose; aap.options.verbose = -1; % ignore error
while isempty(img) && (d < numel(domains))
    d = d + 1;
    domain = domains{d};
    [img, md5, streamdesc] = aas_getfiles_bystream(aap,domain,index{2}(1:end-d+1),streamname,source);
end
aap.options.verbose = verbose0;

if isempty(img)
    aas_log(aap,true,sprintf('%s stream %s not found',source,streamname));
end
end