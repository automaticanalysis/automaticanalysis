% FSL-like overlay (wiht structural) for registration diagnostics
% aas_fsl_overlay(aap,subj[,sess],'input'[,'input'],reslice)
% index must be lower steamlevel
% input can be:
%   - stream
%   - path to image
% 2nd image will be resliced if needed (put lower resolution image first)

function aas_fsl_overlay(aap,varargin)

% parse arguments
index{1} = varargin{1};
if isnumeric(varargin{2})
    index{2} = varargin{2};
    iarg0 = 4;
else
    iarg0 = 3;
end
if isnumeric(index{1})
    switch numel(index)
        case 1 % subject level
            index{2} = index{1};
            index{1} = 'subject';
        case 2 % session level
            index{2} = cell2mat(index);
            index{1} = aap.tasklist.currenttask.domain;
    end
end

for iarg = iarg0:nargin
    if ~ischar(varargin{iarg-1}), break; end
    image{iarg-iarg0+1} = varargin{iarg-1};
    if ~isempty(strfind(image{iarg-iarg0+1},'.')) % path to image
        imagename{iarg-iarg0+1} = basename(image{iarg-iarg0+1});
    else % stream
        imagename{iarg-iarg0+1} = image{iarg-iarg0+1};
        img = '';
        email0 = aap.options.email; aap.options.email = ''; % silence;
        try img = aas_getfiles_bystream(aap,index{:},image{iarg-iarg0+1},'output'); catch, end
        if isempty(img)
            try img = aas_getfiles_bystream(aap,'subject',index{2}(1),image{iarg-iarg0+1},'output'); catch, end
        end
        if isempty(img)
            try img = aas_getfiles_bystream(aap,'study',[],image{iarg-iarg0+1},'output'); catch, end
        end
        aap.options.email = email0;
        if ~isempty(img)
            image{iarg-iarg0+1} = img;
            aas_log(aap,0,'Ignore previous error message(s)!');
        else
            aas_log(aap,1,sprintf('stream %s not found',imagename{iarg-iarg0+1}));
        end
    end
    image{iarg-iarg0+1} = [image{iarg-iarg0+1} ',1']; % ensure single volume
end

% Obtain images
flags = aap.spm.defaults.coreg.write; flags.which = [2 0];
spm_reslice(image,flags)
image{1} = fullfile(fileparts(image{1}),[aap.spm.defaults.coreg.write.prefix basename(image{1}) '.nii']);
image{2} = fullfile(fileparts(image{2}),[aap.spm.defaults.coreg.write.prefix basename(image{2}) '.nii']);

% set path
subj_dir = aas_getpath_bydomain(aap,index{:});

% Create FSL-like overview
iP = fullfile(subj_dir,['diagnostic_' aap.tasklist.main.module(aap.tasklist.currenttask.modulenumber).name '_' imagename{1} '2' imagename{2}]);
aas_runfslcommand(aap,sprintf('slices %s %s -s 3 -o %s.gif',image{2},image{1},iP));
[img,map] = imread([iP '.gif']); s3 = size(img,1)/3;
img = horzcat(img(1:s3,:,:),img(s3+1:2*s3,:,:),img(s3*2+1:end,:,:));
imwrite(img,map,[iP '.jpg']); delete([iP '.gif']);

iP = fullfile(subj_dir,['diagnostic_' aap.tasklist.main.module(aap.tasklist.currenttask.modulenumber).name '_' imagename{2} '2' imagename{1}]);
aas_runfslcommand(aap,sprintf('slices %s %s -s 3 -o %s.gif',image{1},image{2},iP));
[img,map] = imread([iP '.gif']); s3 = size(img,1)/3;
img = horzcat(img(1:s3,:,:),img(s3+1:2*s3,:,:),img(s3*2+1:end,:,:));
imwrite(img,map,[iP '.jpg']); delete([iP '.gif']);
% Clean
delete(image{1}); delete(image{2});
end