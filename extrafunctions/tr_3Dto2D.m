function t = tr_3Dto2D(varargin) %a,nstart,nend,nlag,mx,my
switch nargin
    case 1
        a = varargin{1};
        nstart = 1;
        nend = size(a,3);
        nlag = 1;
        ns = nend - nstart + 1;
        mx = sqrtb(ns/nlag);
        my = mx;
    case 3
        a = varargin{1};
        nstart = 1;
        nend = size(a,3);
        nlag = 1;
        ns = nend - nstart + 1;
        mx = varargin{2};
        my = varargin{3};
    case 4
        a = varargin{1};
        nstart = varargin{2};
        nend = varargin{3};
        nlag = varargin{4};
        if ~nstart 
            nstart = 1;
            nend = size(a,3);        
        end
        ns = nend - nstart + 1;
        mx = sqrtb(ns/nlag);
        my = mx;
    case 5
        a = varargin{1};
        nstart = varargin{2};
        nend = varargin{3};
        nlag = varargin{4};
        if ~nstart 
            nstart = 1;
            nend = size(a,3);        
        end
        ns = nend - nstart + 1;
        mx = varargin{5};
        my = ceil(ns/nlag/mx);
    otherwise
        a = varargin{1};
        nstart = varargin{2};
        nend = varargin{3};
        nlag = varargin{4};
        ns = nend - nstart + 1;
        mx = varargin{5};
        my = varargin{6};
end
fovy = size(a,1);
fovx = size(a,2);
ise = false;
n = nstart-nlag;
t = zeros(my*fovy, mx*fovx);
for y = 1:my
    for x = 1:mx
        n = n + nlag;
        if n > nend
            ise = true;
            break;
        end
        t((y-1)*fovy+1:y*fovy,(x-1)*fovx+1:x*fovx) = a(:,:,n);
    end
    if ise
        break;
    end
end
end