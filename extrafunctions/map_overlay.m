% Transparent overlay of an SPM on a background (bg) image.
% Positive and negative values presented with different colormaps.
% Tibor Auer MRC CBU Cambridge 2012-2013

function [img0 cmap v] = map_overlay(bg,stat,trans)
% Gray - Blue - Yellow-Red (Stronger change is darker)
% cmap = vertcat(gray(128), create_grad([0 0 1],[1 1 1],64),...
    % create_grad([1 1 1],[1 1 0],32),create_grad([1 1 0],[1 0 0],32));

% Gray - Blue - Red-Yellow (Stronger change is brighter)
cmap = vertcat(gray(128), create_grad([1 1 1],[0 0 1],64),...
    create_grad([1 0 0],[1 1 0],32),create_grad([1 1 0],[1 1 1],32));

if nargin < 3, trans = 1; end

if any(bg(:))
    fa = stat.*(stat>0);
    fd = stat.*(stat<0);
    f = zeros(size(bg,1),size(bg,2),size(bg,3));
    if numel(fa(fa>0))
        [fa va] = histc2D(fa,194,256);
        v(194:256) = va;
        f = fa;
    end
    if numel(fd(fd<0))
        [fd vd] = histc2D(fd,130,191);
        v(130:191) = vd;
        f = f + fd;
    end
    img = histc2D(bg,1,127);
    
    for z = 1:size(img,3)
        f0(:,:,z,:) = ind2rgb(f(:,:,z),cmap);
        img0(:,:,z,:) = ind2rgb(img(:,:,z),cmap);
    end
    
    for x = 1:size(img,1)
        for y = 1:size(img,2)
            for z = 1:size(img,3)
                if ~isnan(f(x,y,z)) && f(x,y,z)
                    cr = create_grad(img0(x,y,z,:),f0(x,y,z,:),11);
                    cr(cr<0) = -cr(cr<0); cr(cr>1) = 1;
                    img0(x,y,z,:) = cr(round(trans*10)+1,:);
                end
            end
        end
    end
end
v = v';