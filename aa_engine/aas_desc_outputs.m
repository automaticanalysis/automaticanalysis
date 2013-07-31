% aa describe files that form an output stream 
% Preferred new syntax:
%  function [aap]=aas_desc_outputs(aap,domain,indices,streamname,outputs)
%  e.g., 
%   aas_desc_outputs(aap,'subject',{1},'epi',fns)
%   aas_desc_outputs(aap,'session',{1,1},'epi',fns)
%   aas_desc_outputs(aap,'searchlight',{4,2,77},'epi',fns)
%     [subject 4, session 2, searchligh 77]
%
% Old syntax (still alllowed): 
%  function [aap]=aas_desc_outputs(aap,[subject,[session]],streamname,outputs)
% Outputs may either be specified using a relative or absolute path. If the
% latter, the path is trimmed to make it relative.
% They may be provided in a cell array or as a matrix

function [aap]=aas_desc_outputs(aap,varargin)
global aaworker
% Output stream descriptor
osd=[];

streamname=varargin{end-1};
pos=find(streamname=='.');
if (~isempty(pos))
    streamname=streamname(pos(end)+1:end);
end;
outputs=varargin{end};
streamnme=sprintf('stream_%s_outputfrom_%s.txt',streamname,aap.tasklist.currenttask.name);

osd.name=streamname;
osd.filename=streamnme;

streamdesc='';

switch(nargin)
    case 3
        localroot=aap.acq_details.root;
    case 4
        i=varargin{1};
        localroot=aas_getsubjpath(aap,i);
        [pth nme ext]=fileparts(localroot);
        streamdesc=sprintf(' %s ',[nme ext]);
    case 5
        if ischar(varargin{1})
            localroot=aas_getpath_bydomain(aap,varargin{1},varargin{2});
            streamdesc=localroot;
        else
            i=varargin{1};
            j=varargin{2};
            localroot=aas_getsesspath(aap,i,j);
            [pth nme1 ext1]=fileparts(localroot);
            [pth nme2 ext2]=fileparts(pth);
            streamdesc=sprintf(' %s %s ',[nme2 ext2],[nme1 ext1]);
        end;
end;

osd.desc=streamdesc;

% If outputs provided as cell array, reformat
if (iscell(outputs))
    outputs_asmatrix=[];
    for d=1:length(outputs)
        outputs_asmatrix=strvcat(outputs_asmatrix,outputs{d});
    end;
    outputs=outputs_asmatrix;
end;

osd.numoutputs=size(outputs,1);


descriptor=fullfile(localroot,streamnme);


% Trim absolute path if it has been provided
trimmedoutputs={};
for d=1:size(outputs,1)
    fle=outputs(d,:);
    if (length(fle)>length(localroot) && strcmp(fle(1:length(localroot)),localroot))
        fle=fle(length(localroot)+1:end);
        if (fle(1)==filesep)
            fle=fle(2:end);
        end;
    end;
    trimmedoutputs{d}=deblank(fle);
end;

% Calculate MD5
[aap md5_base64]=aas_md5(aap,trimmedoutputs,localroot);

% Write stream descriptor
fid=fopen(descriptor,'w');
fprintf(fid,'MD5\t%s\n',md5_base64);
for d=1:length(trimmedoutputs)
    fprintf(fid,'%s\n',trimmedoutputs{d});
end;
fclose(fid);


% And propgate to remote filesystem?

switch(aap.directory_conventions.remotefilesystem)
    case 's3'
        aas_log(aap,false,sprintf('Copying %d files to remote filesystem',size(outputs,1)));
        switch(nargin)
            case 3
                s3root=aas_getstudypath(aap,'s3');
            case 4
                i=varargin{1};
                s3root=aas_getsubjpath(aap,i,'s3');
            case 5
                i=varargin{1};
                j=varargin{2};
                s3root=aas_getsesspath(aap,i,j,'s3');
        end;
        s3_copyto_filelist(aap,localroot,streamnme,aaworker.bucket,s3root);
        s3_copyto_filelist(aap,localroot,trimmedoutputs,aaworker.bucket,s3root);
        attr=[];
        % Now write to [user]_streams domain
        attr.numfiles=size(outputs,1);
        attr.md5_base64=md5_base64;
        
        streamentry=[aaworker.bucket '|' s3root '/' streamnme];
        sdb_put_attributes(aap,aaworker.streamtablename,streamentry,attr);
        
        fprintf('Written sdb to %s called %s\n',aaworker.streamtablename,streamentry);

        osd.s3root=s3root;
        aaworker.outputstreams=[aaworker.outputstreams osd];
end;


aas_log(aap,false,sprintf(' output stream %s %s written with %d file(s)',streamname,streamdesc,size(outputs,1)),aap.gui_controls.colours.outputstreams);
