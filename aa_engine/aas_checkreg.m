% FSL-like overlay + SPM orthview video (wiht structural) for registration diagnostics
% aas_checkreg(aap,subj[,sess],'input'[,'input'],reslice)
% index must be lower steamlevel
% input can be:
%   - stream
%   - path to image
% 2nd image will be resliced if needed (put lower resolution image first)

function aas_checkreg(aap,varargin)

%% Parse arguments
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
            if ~isempty(strfind(aap.tasklist.currenttask.domain,'session'))
                index{1} = aap.tasklist.currenttask.domain;
            else
                index{1} = 'session';
            end
    end
end

for iarg = iarg0:nargin
    if ~ischar(varargin{iarg-1}), break; end
    image{iarg-iarg0+1} = varargin{iarg-1};
    if ~isempty(strfind(image{iarg-iarg0+1},'.')) % path to image
        imagename{iarg-iarg0+1} = strtok_ptrn(basename(image{iarg-iarg0+1}),'-0');
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
    image{iarg-iarg0+1} = strcat(image{iarg-iarg0+1},repmat(',1',[size(image{iarg-iarg0+1},1) 1])); % ensure single volume
end

%% SPM
if ~isfield(aap.tasklist.currenttask.settings,'diagnostic') ||...
        isstruct(aap.tasklist.currenttask.settings.diagnostic) ||...
        aap.tasklist.currenttask.settings.diagnostic
    % This will only work for 1-7 overlays
    OVERcolours = {[1 0 0], [0 1 0], [0 0 1], ...
        [1 1 0], [1 0 1], [0 1 1], [1 1 1]};

    % One-by-one
    for i = 1:size(image{1},1)
%         spm_check_registration(image{2});
%         spm_orthviews('addcolouredimage',1,image{1}(i,:), OVERcolours{i})
        spm_check_registration(char(image{2},image{1}(i,:)));
        % switch on Contours
        global st;
        for v = 1:2
            try
                [h f] = get_contextmenu_cb(st.vols{v}.ax{1}.ax,'Contour|Display|all but');
                f(h,[]);
                hM = get_contextmenu_cb(st.vols{v}.ax{1}.ax,'Contour');
                UDc = get(hM,'UserData'); UDc.nblines = 1; set(hM,'UserData',UDc); % narrow
                spm_ov_contour('redraw',v,{});
            catch
            end;
        end
        aas_checkreg_avi(aap, index, 0, ['_' strtok_ptrn(basename(image{1}(i,:)),'-0')]);
        close(1); clear global st;
    end
    if i > 1
        % Summary
        spm_check_registration(image{2});
        for i = 1:size(image{1},1)
            spm_orthviews('addcolouredimage',1,image{1}(i,:), OVERcolours{i})
        end
        aas_checkreg_avi(aap, index, 0,'_summary');
    end
end
% %% FSL
% % Reslice images
% flags = aap.spm.defaults.coreg.write; flags.which = [2 0];
% spm_reslice(image,flags)
% for i = 1:numel(image)
%     for j = 1:size(image{i},1)
%         rimage{i}{j} = fullfile(fileparts(image{i}(j,:)),[aap.spm.defaults.coreg.write.prefix basename(image{i}(j,:)) '.nii']);
%     end
% end
% % binarise if specified
% if isfield(aap.tasklist.currenttask.settings,'PVE') && ~isempty(aap.tasklist.currenttask.settings.PVE)
%     for i = 1:numel(rimage{1})
%         inf = spm_vol(deblank(rimage{1}{i}));
%         Y = spm_read_vols(inf);
%         Y = Y>=aap.tasklist.currenttask.settings.PVE;
%         nifti_write(deblank(rimage{1}{i}),Y,'Binarized',inf)
%     end
% end
% 
% % set path
% subj_dir = aas_getpath_bydomain(aap,index{:});
% 
% % Create FSL-like overview
% for i = 1:numel(rimage{1})
%     image{1} = rimage{1}{i};
%     numel(aap.spm.defaults.coreg.write.prefix)
%     imagename{1} = strtok_ptrn(basename(rimage{1}{i}),'-0');
%     imagename{1} = imagename{1}(1+numel(aap.spm.defaults.coreg.write.prefix):end);
%     image{2} = rimage{2}{1};
%     
%     iP = fullfile(subj_dir,['diagnostic_' aap.tasklist.main.module(aap.tasklist.currenttask.modulenumber).name '_' imagename{1} '2' imagename{2}]);
%     aas_runfslcommand(aap,sprintf('slices %s %s -s 3 -o %s.gif',image{2},image{1},iP));
%     [img,map] = imread([iP '.gif']); s3 = size(img,1)/3;
%     img = horzcat(img(1:s3,:,:),img(s3+1:2*s3,:,:),img(s3*2+1:end,:,:));
%     imwrite(img,map,[iP '.jpg']); delete([iP '.gif']);
%     
%     iP = fullfile(subj_dir,['diagnostic_' aap.tasklist.main.module(aap.tasklist.currenttask.modulenumber).name '_' imagename{2} '2' imagename{1}]);
%     aas_runfslcommand(aap,sprintf('slices %s %s -s 3 -o %s.gif',image{1},image{2},iP));
%     [img,map] = imread([iP '.gif']); s3 = size(img,1)/3;
%     img = horzcat(img(1:s3,:,:),img(s3+1:2*s3,:,:),img(s3*2+1:end,:,:));
%     imwrite(img,map,[iP '.jpg']); delete([iP '.gif']);
%     % Clean
%     delete(image{1});
% end
% delete(image{2});
end