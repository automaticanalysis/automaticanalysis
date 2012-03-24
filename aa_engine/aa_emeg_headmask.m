function aa_emeg_headmask(varargin)
% aa_emeg_headmask([filename [,'display']])
% filename - full path to T1 structural (prompts if not given or not valid)
% 'display' - option to show original and masked volumes 
%
% Attempts to remove vitamin E capsules and other small features outside the head.
% Saves mask and masked version of volume with suffix '_mask' and '_masked'.
%
% Danny Mitchell 07/02/08

if ~isempty(varargin) && exist(varargin{1},'file')
    try V=spm_vol(varargin{1});
    catch delete(varargin{1});
        error('Failed to load %s; file may be corrupt?',varargin{1})
    end
else
    [V sts]=spm_select(1,'image','Select T1 weighted MRI','',pwd,'.*');
    if sts; V=spm_vol(V);
    else return
    end
end

%%%%%%%%%%%% do it
tic
fprintf('\nThresholding...')
% 1) binarise the image
Y=spm_read_vols(V);
BW=Y>mean(Y(:));

% fprintf('\nEroding...')
% 2) this might help seperate capsules very close to the scalp?
% BW=imerode(BW,ones(2,2,2));

fprintf('\b\b\b\b\b\bed. Opening...')
% 3) remove small objects
% CBU vitamin capsules ~1100mm^3
P=spm_imatrix(V.mat);
scal=abs(prod(P(7:9)));
BW=bwareaopen(BW,floor(1150*scal)); % 1200 sometimes removed an eyeball

fprintf('\b\b\b\b\b\bed. Dilating...')
% BW=spm_dilate(BW); crashes - not sure why
% 4) try to close the oesophagus etc
BW=imdilate(BW,ones(6,6,6));

fprintf('\b\b\b\b\b\bed. Filling...')
% 5) fill any holes in the head
warning off all
BW=imfill(BW,'holes');
warning on all

fprintf('\b\b\b\b\b\bed. Eroding...')
% 6) largely undo effect of 4, but avoid eroding scalp
BW=imerode(BW,ones(4,4,4));

fprintf('\b\b\b\b\b\bed. Masking...')
M=BW.*Y;
fname=V.fname;
V.fname=[fname(1:end-4) '_masked.nii'];
spm_write_vol(V,M);
mdfname=V.fname;
V.fname=[fname(1:end-4) '_mask.nii'];
spm_write_vol(V,double(BW));
mfname=V.fname;
clear BW

fprintf('\b\b\b\b\b\bed (%g seconds).',toc)

if nargin~=1
    fprintf(' Displaying...')
    diff=Y-M;
    [m ind]=max(diff(:));
    [x y z]=ind2sub(size(Y),ind);
    xyzmm=V.mat*[x y z 1]';
    V.fname=strrep(mfname,'mask','maskdiff');
    spm_write_vol(V,diff);
    spm_check_registration(char(fname,mdfname,mfname,V.fname));
    spm_orthviews('reposition',xyzmm(1:3));
    ann=annotation('textbox',[0 .94 1 .06],'String',{sprintf('File: %s.',fname),'Raw T1 structural, masked output, difference, and mask.',sprintf('Date:...%s',datestr(now,0))},'edge','none');
    try % print to figures directory if it exists (aa directory format)
        print('-djpeg90',regexprep(mdfname,{'structurals','masked.nii'},{'figures','masking.jpg'}));
    catch % otherwise print to same directory as structural
        print('-djpeg90',regexprep(mdfname,'masked.nii','masking.jpg'));
    end
    delete(ann)
    fprintf('\b\b\b\b\b\bed. \n')
end

return