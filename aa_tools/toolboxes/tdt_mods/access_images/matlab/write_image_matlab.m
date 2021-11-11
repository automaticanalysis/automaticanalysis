function varargout = write_image_matlab(varargin)

%   write_image:
%       inputs: header, volume (X x Y x Z)
%       output: written header (normally not needed)

[hdr, img] = deal(varargin{:});
if ~isequal(hdr.dim(1:ndims(img)),size(img)), error('size of the volume and in header do nor match'); end
fname = strsplit(hdr.fname,','); 
if numel(fname) < 2, fname{2} = 'result'; end
if numel(fname) < 3, fname{3} = 1; end

if fname{3} > 1
    dat = load(fname{1},fname{2});
    dat.(fname{2})(:,:,:,fname{3}) = img;
    img = dat.(fname{2});
end
createVar(fname{2},img);
if ~exist(fname{1},'file')
    save(fname{1},fname{2});
else
    save(fname{1},fname{2},'-append');
end
varargout{1} = hdr;
end

function createVar(name,val)
assignin('caller',name,val);
end