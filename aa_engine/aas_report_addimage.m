% aa report helper function - adds an image to the HTML report
% Rhodri Cusack MRC CBU Cambridge Mar 2006
% Now it uses aas_report_add to support multi-level reporting
% Tibor Auer MRC CBU Cambridge 2012-2013

function aap=aas_report_addimage(varargin)

aap = varargin{1};
for iargin = 1:nargin
    if ischar(varargin{iargin}) && ~isempty(findstr(varargin{iargin},'jpg')), break; end
end
imgpath = varargin{iargin};

if nargin > iargin
    scaling=varargin{nargin};
else
    scaling=1;    
end;

if (~strcmp(imgpath(1:length(aap.acq_details.root)),aap.acq_details.root))
    aas_log(aap,1,sprintf('Cannot relate file %s to directory root %s\n'));
end;

if (scaling==1)
    varargin{iargin} = ['<img src="' imgpath '">'];    
else
    tmp=imread(imgpath);
    varargin{iargin} = ['<img src="' imgpath '" ' sprintf('width=%d height=%d',round(size(tmp,2)*scaling),round(size(tmp,1)*scaling)) '>'];
end;
aap = aas_report_add(varargin{1:iargin});
aap.report.numdependencies=aap.report.numdependencies+1;
aap.report.dependency{aap.report.numdependencies}=imgpath;

