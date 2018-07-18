function DKIadmetrics_vu(DKIcomps_nii,sv_mask)

archstr = computer('arch');
st_comp=archstr(1:3);
if strcmp(st_comp,'win')
    sc='\';
else
    sc='/';
end

% Abrir dados
Vol=load_untouch_nii(sv_mask);
Mask=double(Vol.img);
Vol=load_untouch_nii(DKIcomps_nii);
DKIcomps=double(Vol.img);
[filepth, nme, daa]=fileparts(DKIcomps_nii);

%inicializacao de variaveis
[N1,N2,N3, daa]=size(DKIcomps);
ZER1=zeros(N1,N2,N3);

AWF=ZER1;
Da=zeros(N1,N2,N3,6);
De=Da;
tortu=ZER1;
DAa=ZER1;
DAe=ZER1;
DRe=ZER1;
Kmax=ZER1;

load('Dir125.mat')
for k=1:N3
    for j=1:N2
        for i=1:N1
            if(Mask(i,j,k)==1)
                D=[DKIcomps(i,j,k,1) DKIcomps(i,j,k,2) DKIcomps(i,j,k,3) ...
                    DKIcomps(i,j,k,4) DKIcomps(i,j,k,5) DKIcomps(i,j,k,6)];
                W=[DKIcomps(i,j,k,7) DKIcomps(i,j,k,8) DKIcomps(i,j,k,9) DKIcomps(i,j,k,10) DKIcomps(i,j,k,11) ...
                    DKIcomps(i,j,k,12) DKIcomps(i,j,k,13) DKIcomps(i,j,k,14) DKIcomps(i,j,k,15) DKIcomps(i,j,k,16) ...
                    DKIcomps(i,j,k,17) DKIcomps(i,j,k,18) DKIcomps(i,j,k,19) DKIcomps(i,j,k,20) DKIcomps(i,j,k,21)];
                [AWFi,Dai,Dei,tortui,ADai,ADei,RDei,Kmaxi]=single_population_metrics_ElsF(D,W,V,Nvis,Uvis);
                AWF(i,j,k)=AWFi;
                Da(i,j,k,:)=Dai;
                De(i,j,k,:)=Dei;
                tortu(i,j,k)=tortui;
                DAa(i,j,k)=ADai;
                DAe(i,j,k)=ADei;
                DRe(i,j,k)=RDei;
                Kmax(i,j,k)=Kmaxi;
            end
        end
    end
    disp(k)
end

filepth=[filepth,sc, 'DKIadmetrics',sc];
mkdir(filepth)

AWFnii=Vol;
AWFnii.hdr.dime.dim(5)=1;
AWFnii.hdr.dime.pixdim(5)=0;
AWFnii.hdr.dime.xyz_t=0;
AWFnii.hdr.dime.datatype=16;
AWFnii.img=AWF;
save_untouch_nii(AWFnii,[filepth,sc,nme,'_AWF.nii'])

Danii=Vol;
Danii.hdr.dime.dim(5)=3;
Danii.hdr.dime.pixdim(5)=0;
Danii.hdr.dime.xyz_t=0;
Danii.hdr.dime.datatype=16;
Danii.img=Da;
save_untouch_nii(Danii,[filepth,sc,nme,'_Da.nii'])

Denii=Vol;
Denii.hdr.dime.dim(5)=3;
Denii.hdr.dime.pixdim(5)=0;
Denii.hdr.dime.xyz_t=0;
Denii.hdr.dime.datatype=16;
Denii.img=De;
save_untouch_nii(Denii,[filepth,sc,nme,'_De.nii'])

tortunii=Vol;
tortunii.hdr.dime.dim(5)=1;
tortunii.hdr.dime.pixdim(5)=0;
tortunii.hdr.dime.xyz_t=0;
tortunii.hdr.dime.datatype=16;
tortunii.img=tortu;
save_untouch_nii(tortunii,[filepth,sc,nme,'_tortu.nii'])

DAanii=Vol;
DAanii.hdr.dime.dim(5)=1;
DAanii.hdr.dime.pixdim(5)=0;
DAanii.hdr.dime.xyz_t=0;
DAanii.hdr.dime.datatype=16;
DAanii.img=DAa;
save_untouch_nii(DAanii,[filepth,sc,nme,'_DAa.nii'])

DAenii=Vol;
DAenii.hdr.dime.dim(5)=1;
DAenii.hdr.dime.pixdim(5)=0;
DAenii.hdr.dime.xyz_t=0;
DAenii.hdr.dime.datatype=16;
DAenii.img=DAe;
save_untouch_nii(DAenii,[filepth,sc,nme,'_DAe.nii'])

DRenii=Vol;
DRenii.hdr.dime.dim(5)=1;
DRenii.hdr.dime.pixdim(5)=0;
DRenii.hdr.dime.xyz_t=0;
DRenii.hdr.dime.datatype=16;
DRenii.img=DRe;
save_untouch_nii(DRenii,[filepth,sc,nme,'_DRe.nii'])

Kmaxnii=Vol;
Kmaxnii.hdr.dime.dim(5)=1;
Kmaxnii.hdr.dime.pixdim(5)=0;
Kmaxnii.hdr.dime.xyz_t=0;
Kmaxnii.hdr.dime.datatype=16;
Kmaxnii.img=Kmax;
save_untouch_nii(Kmaxnii,[filepth,sc,nme,'_Kmax.nii'])


