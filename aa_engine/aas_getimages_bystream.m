% Automatic Analysis function to retrieve image lists
% function [imagefns md5]=aas_getimages_bystream(aap,i,j,streamname,varargin)
%  select images from subject i session j
%
% aas_getimages_bystream is intended for EPI images. For other file types, see
% aas_getfiles_bystream.
%
% Rhodri Cusack MRC CBU Cambridge, Nov 2005-2012
% Tibor Auer MRC CBU Cambridge 2012-2013

function [imagefns md5]=aas_getimages_bystream(aap,i,j,streamname,varargin)

% where to look for the data - in the input roots, or the output root
useoutput=true;

v=varargin;
filterimagenumbers=false;

streamname=aas_remapstreamname(aap,streamname,true);

if (~isempty(v))
    if (ischar(v{1}))
        if (strcmp(v{1},'inputroots'))
            useoutput=false;
        end;
        v=v(2:end);
    end;
    if (~isempty(v))
        minimg=v{1};
        v=v(2:end);
        filterimagenumbers=true;
    end;
    if (~isempty(v))
        maximg=v{1};
        filterimagenumbers=true;
    end;
end;

if (~exist('minimg','var'))
    minimg=0;
    maximg=inf;
else
    if (~exist('maximg','var')), maximg=minimg; end;
end;

imagefns=[];

if (useoutput)
    roots={aap.acq_details.root};
else
    roots=aap.acq_details.inputroots;
end;

for inputrootnumber=1:length(roots)
    if (useoutput)
        sesspath=aas_getsesspath(aap,i,j);
    else
        sesspath=aas_getsesspath(aap,i,j,inputrootnumber);
    end;
    
    % Load in all images by stream name
    allimages=[];
    inpstreamdesc=fullfile(sesspath,sprintf('stream_%s_inputto_%s.txt',streamname,aap.tasklist.currenttask.name));
    
    if (~exist(inpstreamdesc,'file'))
        % look for fully qualified inputs
        fqi_filter=fullfile(sesspath,sprintf('stream_*.%s_inputto_%s.txt',streamname,aap.tasklist.currenttask.name));
        fqi=dir(fqi_filter);
        if (length(fqi)>1)
            aas_log(aap, true, sprintf('Found more than one stream matching filter %s - try fully qualifying stream inputs in this module?',fqi_filter));
        elseif (length(fqi)==1)
            inpstreamdesc=fullfile(sesspath,fqi(1).name);
        end;
    end;
    
    fid=fopen(inpstreamdesc,'r');
	% Give error if not found [TA]
    if fid == -1
        aas_log(aap, true, sprintf('%s not found!',inpstreamdesc));
    end
    
    ind=0;
    
    % There should be an MD5 at the top
    lne=fgetl(fid);
    if ((length(lne)>3) && strcmp(lne(1:3),'MD5'))
        md5=lne;
    else
        ind=ind+1;
        allimages(ind).name=lne;
    end;
    
    % Now read in the files
    while (~feof(fid))
        ind=ind+1;
        allimages(ind).name=fgetl(fid);
    end;
    
    fclose(fid);
    
    
    pth=sesspath;
    fns=[];
    for d=1:length(allimages)
        [t r]=strtok(allimages(d).name,'-'); % start of series num
        [t r]=strtok(r,'-'); % start of imgnum
        t=strtok(r,'-');
        imgnum=str2num(t);
        if (~filterimagenumbers || (imgnum>=minimg & imgnum<=maximg))
            imagefns=strvcat(imagefns,fullfile(sesspath,allimages(d).name));
        end;
    end;
    
end;