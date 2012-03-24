% aa report helper function - adds an image to the HTML report
% Rhodri Cusack MRC CBU Cambridge Mar 2006
%

function [aap resp]=aas_report_addimage(aap,imgpath,scaling)

if (~exist('scaling','var'))
    scaling=1;
end;

if (~strcmp(imgpath(1:length(aap.acq_details.root)),aap.acq_details.root))
    aas_log(aap,1,sprintf('Cannot relate file %s to directory root %s\n'));
end;

if (aap.acq_details.root(length(aap.acq_details.root))==filesep)
    imgpath=imgpath((length(aap.acq_details.root)+1):length(imgpath));
else
    imgpath=imgpath((length(aap.acq_details.root)+2):length(imgpath));
end;

if (scaling==1)
    aap.report.html=strcat(aap.report.html,['<img src="' imgpath '">']);
else
    tmp=imread(imgpath);
    aap.report.html=strcat(aap.report.html,['<img src="' imgpath '" ' sprintf('width=%d height=%d',round(size(tmp,2)*scaling),round(size(tmp,1)*scaling)) '>']);
end;
aap.report.numdependencies=aap.report.numdependencies+1;
aap.report.dependency{aap.report.numdependencies}=imgpath;

