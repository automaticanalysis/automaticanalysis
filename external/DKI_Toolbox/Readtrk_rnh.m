% function [INFO, TRACTS]=NBL_readtrk_teste(rootname,conf,order)
function [INFO, TRACTS]=Readtrk_rnh(trk)

%======================================
%  reading TRK
%======================================
% 
% tmp=['Batch',num2str(order)];
% 
% trk=[rootname,'_',conf.(sprintf('%s',tmp)).track,'.trk'];

fid=fopen(trk,'r','l');

disp('reading header....')

INFO.id_string     = fread(fid,6,'*char');
INFO.dim           = fread(fid,3,'int16');
INFO.voxel_size    = fread(fid,3,'float32');
INFO.origin        = fread(fid,3,'float32');

INFO.n_scalar      = fread(fid,1,'int16');
% n_scalars is just a tag that tells how many
% scalars (for each track point) are saved 
%in the track file

INFO.scalar_name   = fread(fid,[10 20],'*char');
INFO.n_properties  = fread(fid,1,'int16');
INFO.property_name = fread(fid,[10 20],'*char');
INFO.vox_to_ras    = fread(fid,[4 4],'float32');
% INFO.vox_to_ras =   fread(fid, 16, 'float32');

INFO.reserved      = fread(fid,444,'*char');

INFO.voxel_order   = fread(fid,4,'*char');
INFO.pad2          = fread(fid,4,'*char');
INFO.img_orient_pat= fread(fid,6,'float32');
INFO.pad1          = fread(fid,2,'*char');

INFO.invert_x      = fread(fid,1,'*uchar');
INFO.invert_y      = fread(fid,1,'*uchar');
INFO.invert_z      = fread(fid,1,'*uchar');
INFO.swap_xy       = fread(fid,1,'*uchar');
INFO.swap_yz       = fread(fid,1,'*uchar');
INFO.swap_zx       = fread(fid,1,'*uchar');

INFO.n_count       = fread(fid,1,'int32');
% voxel count is basically the number of voxels in connection with
% the fiber tracks. So there is no multiple counting

INFO.version       = fread(fid,1,'int32');
INFO.hdr_size      = fread(fid,1,'int32');


disp('reading data....')

%INFO.n_count=1102900;

TRACTS=cell(INFO.n_count,1);

% coordinates + #of scalars
datasize=double(3+INFO.n_scalar);

if INFO.n_properties > 0
    disp('WARNING! There are track properties here and are not loaded!')
end


% tic
for t=1:INFO.n_count
    
    %according to the body of the track file
    
   trlen      = double(fread(fid,1,'*int32')); 
   TRACTS{t}  = fread(fid,[datasize trlen],'float32')';
   
          
end
% toc
fclose(fid);

%INFO
end
