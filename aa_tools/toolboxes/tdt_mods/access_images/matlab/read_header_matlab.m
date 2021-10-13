function varargout = read_header_matlab(varargin)

%   read_header: 
%       input: name of volume (full path , 1 x n string)
%       output: header of volume

% fname = strsplit(varargin{1},':');
% [fname, varname] = deal(fname{:});
% if strcmp(fname,varname), varname = ''; end
% 
try
    hdr = struct('fname',{},'dim',{},'mat',{});
    for v = reshape(whos('-file',varargin{1}),1,[])
        if numel(v.size) < 4, v.size(end+1:4) = 1; end
        if numel(v.size) > 4, error('Cannor handle variables with >4 dimensions'); end
        
        for vol = 1:v.size(4)
            hdr(end+1).fname = strjoin({varargin{1} v.name num2str(v.size(4))},',');
            hdr(end).dim = v.size(1:3);
            hdr(end).mat = eye(4);
        end
    end    
    varargout{1} = hdr;
catch %#ok<CTCH>
    disp(lasterr)
    error(['Cannot read header, probably due to incompatibility ',...
        'between image format and analysis software used or ',...
        'because image does not exist.'])
end