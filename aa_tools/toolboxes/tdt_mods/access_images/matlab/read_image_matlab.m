function varargout = read_image_matlab(varargin)

%   read_image:
%       input: header (struct variable generated with read_header)
%       output: image in neurological space (left = left, view from top)

hdr = varargin{1};

if nargin > 1 % specific variable
    selVar = contains({hdr.fname},[',' varargin{2} ',']);
    if ~any(selVar)
        error('Some variable to be selected does not exist in the MAT file. Please check.')
    end
    hdr = hdr(selVar);
end
if nargin > 2 % specific 3D volume        
    if any(varargin{3}>size(hdr))
        error('Some volume index to be selected exceeds the available number of headers provided in the 4D header file. Please check.')
    end
    hdr = hdr(varargin{2});
else
    if length(hdr)>1
        error('More than one header provided without a volume selection index. Unclear which image to pick from.')
    end
end

img = []; countVar = 0;
for h = reshape(hdr,1,[])
    fname = strsplit(h.fname,',');
    dat = load(fname{1},fname{2});
    countVar = countVar + 1;
    img(:,:,:,countVar) = dat.(fname{2})(:,:,:,str2double(fname{3}));
end

varargout{1} = img;