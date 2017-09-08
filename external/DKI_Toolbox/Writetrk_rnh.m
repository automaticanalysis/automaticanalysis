function Writetrk_rnh(INFO,TRACTS,outputTRK)
%======================================
%====== Starting Spaghetti=>TRK =======
%======================================

%FOR NBL/BRC to correct from 1mm interpolation bias.
%NBL_writetrk(TRACTS,INFO,1.5,-1,-1,'test_shift.trk')

%hard coded data parameters
id_string           = ['TRACK'];
dim                 = INFO.dim;
voxel_size          = INFO.voxel_size;
origin              = INFO.origin;

for i=1:20
        scalar_name(i,1:10)=0;
end

% if strcmp(conf.technique,'dti')
%     n_scalars=0;
% elseif strcmp(conf.technique,'sd')
%     
    n_scalars=0;
    
    % scalar names
   % name1=['IdxAFD'];
   % name2=['IdxDWI'];
   % name3=['IdxAng'];

    
%scalar_name(1:6,1)=name1;
%scalar_name(1:6,2)=name2;
%scalar_name(1:6,3)=name3;
   
% end

n_properties=0;
for i=1:10
    property_name(i,1:20)=0;
end

reserved(1:508)     = 0;
voxel_order         = ['LAS'];
pad2                = ['LAS'];

image_orientation   = [0 0 0 0 0 0];
pad1                = 0;

invert_x            = uint8(0);
invert_y            = uint8(0);
invert_z            = uint8(0);
swap_xy             = uint8(0);
swap_yz             = uint8(0);
swap_zx             = uint8(0);

n_count             = length(TRACTS);
version             = INFO.version;
hdr_size            = 1000;

%============================================
%  Cleaning spaghetti (removing double points)
    % and going back to mm space
%============================================

%disp('shifting streamling data...')

tic
for t=1:n_count
    TRACTS{t}(:,1)=(TRACTS{t}(:,1))*voxel_size(1);
    TRACTS{t}(:,2)=(TRACTS{t}(:,2))*voxel_size(2);
    TRACTS{t}(:,3)=(TRACTS{t}(:,3))*voxel_size(3);
end
toc

%======================================
%  Writing TRK
%======================================
fid=fopen(outputTRK,'wb','l');

fwrite(fid,id_string,'char');
fwrite(fid,0,'uint8'); %PAD FLAVIO...
fwrite(fid,dim,'int16');
fwrite(fid,voxel_size,'float32');
fwrite(fid,origin,'float32');

fwrite(fid,n_scalars,'int16');
fwrite(fid,scalar_name,'char');
fwrite(fid,n_properties,'int16');
fwrite(fid,property_name,'uint8');

% fwrite(fid,mtx,'float32');

fwrite(fid,reserved,'uint8');
fwrite(fid,voxel_order,'char');
fwrite(fid,0,'uint8'); %PAD FLAVIO
fwrite(fid,pad2,'char');
fwrite(fid,0,'uint8'); %PAD FLAVIO
fwrite(fid,image_orientation,'float32');
fwrite(fid,pad1,'int16'); %double pad uint8

fwrite(fid,invert_x,'uchar');
fwrite(fid,invert_y,'uchar');
fwrite(fid,invert_z,'uchar');
fwrite(fid,swap_xy,'uchar');
fwrite(fid,swap_yz,'uchar');
fwrite(fid,swap_zx,'uchar');

fwrite(fid,n_count,'int32');
fwrite(fid,version,'int32');
fwrite(fid,hdr_size,'int32');

disp('saving....')

tic
for t=1:n_count
   fwrite(fid,size(TRACTS{t},1),'int32'); 
   fwrite(fid,TRACTS{t}','float32');
end
toc
fclose(fid);
end