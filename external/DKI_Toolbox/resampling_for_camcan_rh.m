function RV=resampling_for_camcan_rh(V,resinitial,resfinal)
% Resamples volume of image V with initial resolution 'resinitial' to a
% resolution 'resfinal'
% Rafael Neto Henriques, 2013

sizInitial=size(V);
sizXinitial=sizInitial(1);
sizYinitial=sizInitial(2);
sizZinitial=sizInitial(3);

positionXvoxel=(1:sizXinitial)*resinitial(1)-resinitial(1)/2;
positionYvoxel=(1:sizYinitial)*resinitial(2)-resinitial(2)/2;
positionZvoxel=(1:sizZinitial)*resinitial(3)-resinitial(3)/2;

FOVx=sizXinitial*resinitial(1);
FOVy=sizYinitial*resinitial(2);
FOVz=sizZinitial*resinitial(3);

Xfinal=resfinal(1)/2:resfinal(1):FOVx;
Yfinal=resfinal(2)/2:resfinal(2):FOVy;
Zfinal=resfinal(3)/2:resfinal(3):FOVz;

if Xfinal(1)<positionXvoxel(1)
    Xfinal(1)=positionXvoxel(1)+eps;
end
if Zfinal(1)<positionZvoxel(1)
    Zfinal(1)=positionZvoxel(1)+eps;
end
if Yfinal(1)<positionYvoxel(1)
    Yfinal(1)=positionYvoxel(1)+eps;
end

if Xfinal(end)>positionXvoxel(end)
    Xfinal(end)=positionXvoxel(end)-eps;
end
if Zfinal(end)>positionZvoxel(end)
    Zfinal(end)=positionZvoxel(end)-eps;
end
if Yfinal(end)>positionYvoxel(end)
    Yfinal(end)=positionYvoxel(end)-eps;
end

[X,Y,Z] = meshgrid(positionXvoxel, positionYvoxel, positionZvoxel);

[Xi,Yi,Zi] = meshgrid(Xfinal, Yfinal, Zfinal);

RV = interp3(X,Y,Z, V,Xi,Yi,Zi,...
    '*cubic');
