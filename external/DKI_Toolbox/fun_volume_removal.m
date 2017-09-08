function [outDWI,outbval,outbvec]=fun_volume_removal(DWI,bval,bvec,num)

archstr = computer('arch');
st_comp=archstr(1:3);
if strcmp(st_comp,'win')
    sc='\';
else
    sc='/';
end

% 1
addpath(['..',sc,'NIFTI_toolbox',sc])

% 2
[pth,source]=fileparts(DWI);
[pth,sbval]=fileparts(bval);
[pth,sbvec]=fileparts(bvec);

sourcepth=[pth,sc];

targetpth=[pth,sc];

source=[source,'.nii'];
sbval=[sbval,'.bval'];
sbvec=[sbvec,'.bvec'];

target=['r_',source];

tbval=['r_',sbval];
tbvec=['r_',sbvec];

% 3
Vs=load_untouch_nii([sourcepth,source]);
idim=Vs.hdr.dime.dim(5);
Vs.hdr.dime.dim(5)=idim-1;
DWImat=Vs.img;
DWImat(:,:,:,num)=[];
Vs.img=DWImat;

% 5
outDWI=[targetpth,target];
save_untouch_nii(Vs,outDWI);

% 6
bl=load([sourcepth,sbval]);
bv=load([sourcepth,sbvec]);
bl(num)=[];
bv(:,num)=[];

% 7
outbval=[targetpth,tbval];
fidl=fopen(outbval,'w');
fprintf(fidl,'%4d ',bl);
fprintf(fidl,'\n');
fclose(fidl);

outbvec=[targetpth,tbvec];
fidv=fopen(outbvec,'w');
fprintf(fidv,'%.14f ',bv(1,:));
fprintf(fidv,'\n');
fprintf(fidv,'%.14f ',bv(2,:));
fprintf(fidv,'\n');
fprintf(fidv,'%.14f ',bv(3,:));
fprintf(fidv,'\n');
fclose(fidv);
