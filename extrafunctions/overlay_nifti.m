function [ errflag,errstring ] = overlay_nifti(nii_fname, template_fname, render_fname, savefile_fname, image_label, brightness)
%
% create a render + three orthogonal overlays for the specified
% nifti(s) using the SPM and aa image functions and save to jpg 
%
% INPUTS
%
%   nii_fname      - cell array of up to 3 nifti to overlay (fullpath if not in working dir)
%   template_fname - structural template (fullpath if not in working dir)
%   render_fname   - render template (fullpath if not in working dir)
%   savefile_fname - output fname (see OUTPUT) (fullpath if not in working dir)
%
% OPTIONAL INPUTS
%
%   image_label    - text used to label figure (empty or not passed = no label)
%   brightness     - render brightness (according to smp_render) but it works more
%                    like contrast: monochrome-ish (0) to foggy (1.0) (default: 0.5)
%
% OUTPUT
%
%   savefile_fname_render.jpg - overlay on rendered brain
%   savefile_fname_01.jpg - overlay on axial slices
%   savefile_fname_02.jpg - overlay on coronal slices
%   savefile_fname_03.jpg - overlay on sagittal slices
%
% EXAMPLE USAGE:
%
%  nii_fname = 'thrT_0001.nii';
%  template_fname = '/Applications/spm12/toolbox/OldNorm/T1.nii';
%  render_fname = '/Applications/spm12/rend/render_single_subj.mat';
%
% NOTES
%
% 1) % If the first char isn't '/', the function will look in under the
% SPM install dir for the template and render nii -- e.g.:
%
%       template_fname = 'toolbox/OldNorm/T1.nii';
%       render_fname = 'rend/render_single_subj.mat';
%
% 2) rendering works for multiple nii but slices don't. Also, the render
% results can be poor if the values in the nii vary substantially (e.g., a
% 0/1 mask and a -5 to 5 tmap). Ostensibly the colors go nii1 = red,
% nii2 = green, nii3 = blue. If more than 3 nii are passed, the function
% barfs and returns.
%
%
% CHANGE HISTORY
%
% 08/2020 [msj] - change default brightness from 0.5 to 0.2
% 10/2019 [MSJ] - modified to handle up to three nii
% 09/2019 [MSJ] - new
%

errflag = 0;
errstring = '';


% ------------------------------------------------------------------------
% sanity checks
% ------------------------------------------------------------------------

if (nargin < 4) 
    errstring = 'Usage: overlay_nifti(nii_fname, template_fname, render_fname, savefile_fname, [image_label, brightness])';
    errflag = 1;
    return;
end
    
if (iscell(nii_fname))
    
    num_nii = numel(nii_fname);

    if (num_nii > 3)
        errstring = 'more than 3 nii passed';
        errflag = 1;
        return;
    end

    for index = 1:num_nii
        p = dir(nii_fname{index});
        if (isempty(p))
            errstring = sprintf('nifti file %s not found', nii_fname{index});
            errflag = 1;
            return;
        end
    end
        
else
    
    num_nii = 1;

    p = dir(nii_fname);
    if (isempty(p))
        errstring = sprintf('nifti file %s not found', nii_fname);
        errflag = 1;
        return;
    end
    
end

p = dir(render_fname);
if isempty(p) && (render_fname(1) ~= '/') 
    render_fname = fullfile(fileparts(which('spm')),render_fname); 
end
     
p = dir(render_fname);
if isempty(p)
    errstring = sprintf('render file %s not found', render_fname);
    errflag = 1;
    return;
end

p = dir(template_fname);
if (template_fname(1) ~= '/') 
    template_fname = fullfile(fileparts(which('spm')),template_fname); 
end

p = dir(template_fname);
if (isempty(p))
    errstring = sprintf('template file %s not found', template_fname);
    errflag = 1;
    return;
end

if (nargin < 5)
    image_label = [];
end


if (nargin < 6)
%     brightness = 0.5;
    % this better for multiple map overlap
    brightness = 0.2;
end


% ------------------------------------------------------------------------
% 1) Render - we can use shortcut in spm_render by just passing
% in the nifti and render image filename
% ------------------------------------------------------------------------

% need to set up global prevrend for some reason...

global prevrend
prevrend = struct('rendfile', render_fname, 'brt',0.5, 'col',eye(3));
out = spm_render(nii_fname, brightness, render_fname);
spm_figure('Close','Graphics');

% squeeze render output into montage and write to file

for i = 1:numel(out), img(1:size(out{i},1),1:size(out{i},2),:,i) = out{i}; end
mon = tr_3Dto2D(squeeze(img(:,:,1,[1 3 5 2 4 6])));
mon(:,:,2) = tr_3Dto2D(squeeze(img(:,:,2,[1 3 5 2 4 6])));
mon(:,:,3) = tr_3Dto2D(squeeze(img(:,:,3,[1 3 5 2 4 6])));
mon = mon(1:size(mon,2)*2/3,:,:);			

if (~isempty(image_label))
    mon = insertInImage(mon, @()text(40,25,strrep(image_label,'_','-')),...
    {'fontweight','bold','color','y','fontsize',16,...
    'linewidth',1,'margin',5,'backgroundcolor','k'});	
end

[ p,n,~ ] = fileparts(savefile_fname);
render_savefile_fname = fullfile(p,[n '_render.jpg']);

imwrite(mon,render_savefile_fname);

  
% ------------------------------------------------------------------------
% 2) three ortho section overlays 
% (this code mostly stolen from aamod_firstlevel_threshold...)
% ------------------------------------------------------------------------

% alas, slices currently can't handle > 1 nii

if (num_nii > 1)
    return;
end

transparency = 0.4;
nth_slice = 3; 

template_header = spm_vol(template_fname);
[ Ytemplate,~ ] = spm_read_vols(template_header);

% work out threshold for template

threshprop=0.10;
ys=sort(Ytemplate(~isnan(Ytemplate)));
bright3=ys(round(length(ys)*0.3));
bright97=ys(round(length(ys)*0.97));
thresh=bright3*(1-threshprop)+bright97*threshprop;
Ytemplate=Ytemplate.*(Ytemplate>thresh);

% need nifti to match template

nifti_header = spm_vol(nii_fname);

if (~isequal(nifti_header.dim, template_header.dim) ||	norm(nifti_header.mat-template_header.mat)>0.01)

    resliceOpts = [];
    resliceOpts.mask = false;
    resliceOpts.mean = false;
    resliceOpts.interp = 1;
    resliceOpts.which = 1;		% DON'T reslice the first image
    resliceOpts.wrap = [0 0 0];	% this is everywhere in aa even though it should be [1 1 0] for MRI
    resliceOpts.prefix = 'r';

    spm_reslice({template_header.fname, nifti_header.fname}, resliceOpts);
    
    [ p,n,e ] = fileparts(nii_fname);
    nifti_header = spm_vol(fullfile(p,[resliceOpts.prefix n e]));

end

[ rYepi,~ ] = spm_read_vols(nifti_header);
        
for a = 0:2 % in 3 axes

    arYepi = shiftdim(rYepi,a);
    aYtemplate = shiftdim(Ytemplate,a);
    
    % Adjust slice selection according to the activation
    
    for iSl = 1:nth_slice 
        iYepi = arYepi(:,:,iSl:nth_slice:end);
        if any(iYepi(:)~=0), break; end
    end
    
    iYepi = img_rot90(iYepi);
    iYtemplate = img_rot90(aYtemplate(:,:,iSl:nth_slice:end));

    [ img,~,~ ] = map_overlay(iYtemplate,iYepi,1-transparency);                                
    mon = tr_3Dto2D(img_tr(img(:,:,:,1),a==2));
    mon(:,:,2) = tr_3Dto2D(img_tr(img(:,:,:,2),a==2));
    mon(:,:,3) = tr_3Dto2D(img_tr(img(:,:,:,3),a==2));

    if (~isempty(image_label))
        mon = insertInImage(mon, @()text(40,25,strrep(image_label,'_','-')),...
        {'fontweight','bold','color','y','fontsize',14,...
        'linewidth',1,'margin',5,'backgroundcolor','k'});	
    end
    
    [ p,n,~ ] = fileparts(savefile_fname);
	slice_fname = fullfile(p,sprintf('%s_%d.jpg',n,a+1));

	imwrite(mon,slice_fname);

end

% delete resliced image if one was created

[ p,n,e ] = fileparts(nii_fname);
try delete(fullfile(p,[resliceOpts.prefix n e])); catch; end
    
end 
 

% ------------------------------------------------------------------------
function fo = img_rot90(fi)
% ------------------------------------------------------------------------

for i = 1:size(fi,3)
    fo(:,:,i) = rot90(fi(:,:,i),1);
end
end


% ------------------------------------------------------------------------
function fo = img_tr(fi,toDo)
% ------------------------------------------------------------------------

if nargin < 2, toDo = true; end
if toDo
    nslice = size(fi,3);
    for i = 1:nslice
        fo(:,:,i) = fliplr(rot90(fi(:,:,i),1));
    end
else
    fo = fi;
end
end
