% Overlay an epi "Yepi" (read with spm_read_vols) onto the template
% Tibor Auer MRC CBU Cambridge 2012-2013

function f = aas_overlay(aap,Yepi,mat)
try nSl = aap.internal.aap_initial.tasksettings.aamod_firstlevel_threshold.overlay.nth_slice; catch, nSl = 3; end;
try tra = aap.internal.aap_initial.tasksettings.aamod_firstlevel_threshold.overlay.transparency; catch, tra = 0.5; end

% Template
Vtemplate=spm_vol(aap.directory_conventions.T1template);
[Ytemplate tXYZ]=spm_read_vols(Vtemplate);
tXYZ=[tXYZ;ones(1,size(tXYZ,2))];
% work out threshold for template
threshprop=0.10;
ys=sort(Ytemplate(~isnan(Ytemplate)));
bright3=ys(round(length(ys)*0.3));
bright97=ys(round(length(ys)*0.97));
thresh=bright3*(1-threshprop)+bright97*threshprop;
Ytemplate=Ytemplate.*(Ytemplate>thresh);

% Resize
iXYZ=mat\tXYZ;
rYepi=spm_sample_vol(Yepi,iXYZ(1,:),iXYZ(2,:),iXYZ(3,:),0);
rYepi=reshape(rYepi,size(Ytemplate));

% Overlay
for iSl = 1:nSl-1 % Adjust slice selection according to the activation
    iYepi = rYepi(:,:,iSl:nSl:end);
    if any(iYepi(:)~=0), break; end
end
iYepi = img_rot90(iYepi);
iYtemplate = img_rot90(Ytemplate(:,:,iSl:nSl:end));
[img cm v] = map_overlay(iYtemplate,iYepi,tra);
mon = tr_RGBtoMontage(img);
f = figure;
montage(mon);
colormap(cm);
cb = colorbar;
yT = [walley(v) find(v<0, 1, 'last' ) find(v>0, 1 ) peak(v)];
set(cb,'YTick',yT,'YTickLabel',v(yT));