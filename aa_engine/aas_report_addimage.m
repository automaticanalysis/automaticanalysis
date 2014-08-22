% aa report helper function - adds an image to the HTML report
% Rhodri Cusack MRC CBU Cambridge Mar 2006
% Now it uses aas_report_add to support multi-level reporting
% Tibor Auer MRC CBU Cambridge 2012-2013

function aap=aas_report_addimage(aap,varargin)

defWIDTH=800;

for iargin = 1:nargin-1
    if ischar(varargin{iargin}) && ~isempty(strfind(varargin{iargin},'jpg')), break; end
end
imgpath = varargin{iargin};

if (nargin-1) > iargin
    scaling=varargin{nargin-1};
else
    scaling=0;    
end;

if (~strcmp(imgpath(1:length(aap.acq_details.root)),aap.acq_details.root))
    aas_log(aap,1,sprintf('Cannot relate file %s to directory root %s\n'),imgpath,aap.acq_details.root);
end;

[p n e] = fileparts(imgpath);
switch e
    case {'.png', '.jpg', '.jpeg'}
        varargin{iargin} = sprintf('<br><cite>%s</cite><br>',basename(imgpath));
        aap = aas_report_add(aap,varargin{1:iargin});
        pfx = '<img src';
        sfx = '><br>';
        tmp=imread(imgpath);
        if ~scaling && (size(tmp,2) <= defWIDTH), scaling = 1; end
    case '.avi'
        pfx = '<a href';
        sfx = sprintf('>Play video: %s</a><br>',n);
        scaling = 1;
    otherwise
        aas_log(aap,1,sprintf('Unknown format: %s\n',e));
end

switch scaling
    case 0
        str = [pfx '="' imgpath '" ' sprintf('width=%d height=%d',defWIDTH,round(defWIDTH/size(tmp,2)*size(tmp,1))) sfx];
    case 1
        str = [pfx '="' imgpath '"' sfx];    
    otherwise
        str = [pfx '="' imgpath '" ' sprintf('width=%d height=%d',round(size(tmp,2)*scaling),round(size(tmp,1)*scaling)) sfx];
end
varargin{iargin} = str;
aap = aas_report_add(aap,varargin{1:iargin});
aap.report.numdependencies=aap.report.numdependencies+1;
aap.report.dependency{aap.report.numdependencies}=imgpath;

