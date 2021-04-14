% FSL-like overlay + SPM orthview video (wiht structural) for registration diagnostics
% aas_checkreg(aap,subj[,sess],'input'[,'input'],reslice)
% index must be lower streamlevel
% input can be:
%   - stream
%   - path to image
% 2nd image will be resliced if needed (put lower resolution image first)
%
% CHANGE HISTORY
%
% 6/17 [MSJ] -- use a unique prefix for reslicing (not 'r')
%


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
    if isempty(strfind(image{iarg-iarg0+1},'.nii')) % stream
        img = aas_getfiles_bystream_multilevel(aap,index{:},image{iarg-iarg0+1},'output');
        if ~isempty(img), image{iarg-iarg0+1} = img; end
    end
    image{iarg-iarg0+1} = strcat(image{iarg-iarg0+1},repmat(',1',[size(image{iarg-iarg0+1},1) 1])); % ensure single volume
end

if ~isfield(aap.tasklist.currenttask.settings,'diagnostic') ||...
        (~isstruct(aap.tasklist.currenttask.settings.diagnostic) && aap.tasklist.currenttask.settings.diagnostic) ||...
        (isstruct(aap.tasklist.currenttask.settings.diagnostic) && aap.tasklist.currenttask.settings.diagnostic.streamind)
    
    % One-by-one
    imgExcl = [];
    for i = 1:size(image{1},1)
        %         spm_check_registration(image{2});
        %         spm_orthviews('addcolouredimage',1,image{1}(i,:), OVERcolours{i})
        isOK = true; try spm_vol(image{1}(i,:)); catch, isOK = false; end
        if ~isOK
            aas_log(aap,false,sprintf('WARNING: trype of file %s is not recognised --> skipping',image{1}(i,:))); 
            imgExcl(end+1) = i;
            continue; 
        end
        spm_check_registration(char(image{2},image{1}(i,:)));
        % switch on Contours
        global st;
        switch spm('ver')
            case {'SPM12b' 'SPM12'}
                for v = 1:2
                    [h f] = get_contextmenu_cb(st.vols{v}.ax{1}.ax,'Contour|Display|all but');
                    f(h,[]);
                    hM = get_contextmenu_cb(st.vols{v}.ax{1}.ax,'Contour');
                    UDc = get(hM,'UserData'); UDc.nblines = 1; set(hM,'UserData',UDc); % narrow
                    spm_ov_contour('redraw',v,{});
                end
        end
        aas_checkreg_avi(aap, index, 0, ['_' strtok_ptrn(basename(image{1}(i,:)),'-0')]);
        close(spm_figure('GetWin','Graphics')); clear global st;
    end
    image{1}(imgExcl,:) = [];
    if i > 1
        % Summary
        OVERcolours = distinguishable_colors(size(image{1},1));

        spm_check_registration(image{2});
        for i = 1:size(image{1},1)
            spm_orthviews('addcolouredimage',1,image{1}(i,:), OVERcolours(i,:))
        end
        aas_checkreg_avi(aap, index, 0,'_summary');
        close(spm_figure('GetWin','Graphics'));
    end
end

end
