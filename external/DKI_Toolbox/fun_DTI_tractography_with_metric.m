function fun_DTI_tractography_with_metric(V1_name,S_name,sT,T_name,Tm_t,...
    tangle,tmin,tmax,NVDSeeds, M_name)
%% Performs DTI based streamline tracking
% Inputs
% 1) V1_name - directory and folder name (with extention) where principal 
% direction of diffusion tensor is saved
% 2) S_name - directory of the map used to select the voxels to seed, 
% i.e. the voxels where tracts will started
% 3) sT - All voxels larger than this threshold on S_name willl be seeded
% for binary 
% 4) FA_name - directory and folder name (with extention) where FA maps is 
% saved. This map this used for the stopping criteria.
% 5) tFA - Value of FA threshold for FA stopping criteria. Tracking will 
% stop if tracts reach a voxel with FA is lower than this threshold.
% 6) tangle - Angle Deviation threshold: tracking stops if new direction 
% have a angle deviation larger than the setted value.
% 7) tmin (cm) - tracts with length lower than this value are excluded
% 8) tmax (cm) - tracts with length higher than this value are excluded
% 9) NVDSeeds - Parameter to set the number of seed per voxel dimention.
% Seeds are distributed on a grid of NVDSeeds3
%
% first version: 04/07/2013
% Rafael Neto Henriques 
%
% Updates:
% 03/12/2014 RNH -> Ajusting stopping criteria for FA: 
%                   point that falls inside FA_name < tFA is also saved
%                   this insure that streamlines that stard in white matter
%                   touches the cortical regions
% 31/01/2015 RNH -> introducing evenly seeding option:
%                   Now the algorithm allow to produce multiple seeds 
%                   inside each voxel. This option is for evenly seeding, 
%                   and NxNxN grid of seeds is generated inside of each
%                   voxel. 
%                   NVDSeeds = N (number of seeds per voxel dimention)
% 11/04/2015 RNH -> Addiding Report metrics, removing code breaks
%%

% Windows, Linux or Mac
archstr = computer('arch');
st_comp=archstr(1:3);
if strcmp(st_comp,'win')
    sc='\';
else
    sc='/';
end
addpath(['..',sc,'NIFTI_toolbox'])

% Define output name
[pth,outname, daa]=fileparts(V1_name);
mkdir([pth, sc,'Tractography'])
outputTRK=[pth,sc,'Tractography', sc, outname, '_tdti.trk'];

% Load data
Vd=load_untouch_nii(V1_name);
Tm=load_untouch_nii(T_name);
Sm=load_untouch_nii(S_name);
VMetric=load_untouch_nii(M_name);

D=Vd.img;
Tmask=Tm.img;
Smask=Sm.img;
Metric=VMetric.img;

Dimen=size(Vd.img);
voxel_size=Vd.hdr.dime.pixdim(2:4);
pd=mean(voxel_size);

Smask(1,:,:)=0; % Two remove volume limits
Smask(end,:,:)=0;
Smask(:,1,:)=0;
Smask(:,end,:)=0;
Smask(:,:,1)=0;
Smask(:,:,end)=0;

Tmask(1,:,:)=0; % Two remove volume limits
Tmask(end,:,:)=0;
Tmask(:,1,:)=0;
Tmask(:,end,:)=0;
Tmask(:,:,1)=0;
Tmask(:,:,end)=0;

% compute normalized thresholds
agv_t=cos(tangle/180*pi); % converting ag_t to analize angular deviations
                          % directly in cross product between vectors
Tmin_t=tmin/pd; % Minimum fiber length (norm)
Tmax_t=tmax/pd; % Maximun fiber length (norm)
step=0.5;       % Distance of forward fiber in one step (norm)

% Define seeds (mark voxels to be seeded)
sj=Dimen(2);SeedJ=1:sj;
si=Dimen(1);SeedI=1:si;
sk=Dimen(3);SeedK=1:sk;
[J,I,K] = meshgrid(SeedJ,SeedI,SeedK);
AllVox(:,3) = K(:)';
AllVox(:,1) = I(:)';% X have to correspond to index I
AllVox(:,2) = J(:)';
SelectVoxels=AllVox(Smask(:)>sT,:); 
NSVoxels=size(SelectVoxels,1);
Tmask_v=Tmask(:);

% Select Seeds out region of cut
SelectVoxels=SelectVoxels(Tmask_v(Smask(:)>sT)>Tm_t,:); 
ANSVoxels=size(SelectVoxels,1); % N actual Seeds
Ptracts=ANSVoxels*NVDSeeds^3;

% Seed voxels according to number of seeds per voxel
Seed=zeros(Ptracts,3);
sx=1/(2*NVDSeeds):1/(NVDSeeds):1-1/(2*NVDSeeds);
sy=sx;
sz=sx;
[sx,sy,sz] = meshgrid(sx,sy,sz);
sx=sx(:);
sy=sy(:);
sz=sz(:);
for s=1:NVDSeeds^3
    Seed(1+(s-1)*ANSVoxels:s*ANSVoxels,1)=SelectVoxels(:,1)-sx(s);
    Seed(1+(s-1)*ANSVoxels:s*ANSVoxels,2)=SelectVoxels(:,2)-sy(s);
    Seed(1+(s-1)*ANSVoxels:s*ANSVoxels,3)=SelectVoxels(:,3)-sz(s);
end

% Report parameters
T_A_Cut=0;
T_R_Cut=0;
T_AR_Cut=0;
T_M1_Cut=0;
T_M2_Cut=0;
E_A_Cut=0;
E_R_Cut=0;
E_AR_Cut=0;
E_M1_Cut=0;
E_M2_Cut=0;
Others=0;
Ethers=0;

t=0;
texcluded=0;

tic
for q = 1:Ptracts;

    Temp_Tracts = [];
    Temp_TractsM = [];
    
    % Starting point at seed
    i = Seed(q,1);
    j = Seed(q,2);
    k = Seed(q,3);
    
    % Seeds corresponding voxel
    Int_i = ceil(i);
    Int_j = ceil(j);
    Int_k = ceil(k);
    
    % inicial point
    s=1;
    Temp_Tracts(s,:)= [i,j,k];
    Temp_TractsM(s,:)= Metric(Int_i,Int_j,Int_k);
    
    %% Start propogation in the first direction
    s=2;
    
    % First propogation is done out of while loop so we can control the
    % direction of propogation
    
    Vout = [D(Int_i,Int_j,Int_k,1);D(Int_i,Int_j,Int_k,2);D(Int_i,Int_j,Int_k,3)];
    i = i + step*Vout(1);
    j = j + step*Vout(2);
    k = k + step*Vout(3);
    Int_i = ceil(i);
    Int_j = ceil(j);
    Int_k = ceil(k);
    Temp_Tracts(s,:) = [i,j,k];
    Temp_TractsM(s,:) = Metric(Int_i,Int_j,Int_k);
    
    Temp_end1_R=true;
    Temp_end1_A=true;
    Temp_end1_M=true;
    
    while (Temp_end1_M && Temp_end1_R && Temp_end1_A );
        
        if Tmask(Int_i,Int_j,Int_k) > Tm_t
            
            s=s+1;
            
            NVout = [D(Int_i,Int_j,Int_k,1);D(Int_i,Int_j,Int_k,2);D(Int_i,Int_j,Int_k,3)];
            Dotproduct = (Vout'*NVout)/(norm(Vout')*norm(NVout'));
            
            if abs(Dotproduct)<agv_t
                Temp_end1_A=false;
            else
                if Dotproduct < 0;
                    NVout = -NVout;
                end
                
                % NVout
                i = i + step*NVout(1);
                j = j + step*NVout(2);
                k = k + step*NVout(3);
                Int_i = ceil(i);
                Int_j = ceil(j);
                Int_k = ceil(k);
                Vout=NVout;
                
                Temp_Tracts(s,:)= [i,j,k];
                Temp_TractsM(s,:)= Metric(Int_i,Int_j,Int_k);
            end
            % Continue loop if pass all restritions
            if s*step >= Tmax_t;% If a single fiber segment bigger than length threshold, stopping fiber tracking
                Temp_end1_M = false;
            end
        else
            Temp_end1_R = false;
        end
    end
    
    if Temp_end1_M
        %% In the second way
        Temp_Tracts = flipud(Temp_Tracts);
        Temp_TractsM =  flipud(Temp_TractsM);
        
        % Reset starting point of single fiber segment in the second direction
        i = Seed(q,1);
        j = Seed(q,2);
        k = Seed(q,3);
        Int_i = ceil(i);
        Int_j = ceil(j);
        Int_k = ceil(k);
        
        Vout = -[D(Int_i,Int_j,Int_k,1);D(Int_i,Int_j,Int_k,2);D(Int_i,Int_j,Int_k,3)];
        i = i + step*Vout(1);
        j = j + step*Vout(2);
        k = k + step*Vout(3);
        Int_i = ceil(i);
        Int_j = ceil(j);
        Int_k = ceil(k);
        
        Temp_Tracts(s,:)= [i,j,k];
        Temp_TractsM(s,:)= Metric(Int_i,Int_j,Int_k);
        
        Temp_end2_R=true;
        Temp_end2_A=true;
        Temp_end2_M=true;
        
        while (Temp_end2_R && Temp_end2_A && Temp_end2_M);
            
            if Tmask(Int_i,Int_j,Int_k) > Tm_t
                
                s=s+1;
                
                NVout = [D(Int_i,Int_j,Int_k,1);D(Int_i,Int_j,Int_k,2);D(Int_i,Int_j,Int_k,3)];
                Dotproduct = (Vout'*NVout)/(norm(Vout')*norm(NVout'));
                
                if abs(Dotproduct)<agv_t
                    Temp_end2_A=false;
                else
                    if Dotproduct < 0;
                        NVout = -NVout;
                    end
                    
                    % NVout
                    i = i + step*NVout(1);
                    j = j + step*NVout(2);
                    k = k + step*NVout(3);
                    Int_i = ceil(i);
                    Int_j = ceil(j);
                    Int_k = ceil(k);
                    Vout=NVout;
                    
                    Temp_Tracts(s,:)= [i,j,k];
                    Temp_TractsM(s,:)= Metric(Int_i,Int_j,Int_k);
                end
                % Continue loop if pass all restritions
                if s*step >= Tmax_t;% If a single fiber segment bigger than length threshold, stopping fiber tracking
                    Temp_end2_M = false;
                end
            else
                Temp_end2_R = false;
            end
        end
    end
    
    if s*step>Tmin_t
        t=t+1;
            Tracts{t,1} = Temp_Tracts;
            TractsM{t,1} = Temp_TractsM;
        
        % report
        if ~Temp_end1_A && ~Temp_end2_A
            T_A_Cut=T_A_Cut+1;
        elseif ~Temp_end1_A && ~Temp_end2_R
            T_AR_Cut=T_AR_Cut+1;
        elseif ~Temp_end1_R && ~Temp_end2_A
            T_AR_Cut=T_AR_Cut+1;
        elseif ~Temp_end1_R && ~Temp_end2_R
            T_R_Cut=T_R_Cut+1;
        elseif ~Temp_end1_M
            T_M1_Cut=T_M1_Cut+1;
        elseif ~Temp_end2_M
            T_M2_Cut=T_M2_Cut+1;
        else
            Others=Others+1;
        end
        
    else
        texcluded=texcluded+1;
        
        % report
        if ~Temp_end1_A && ~Temp_end2_A
            E_A_Cut=E_A_Cut+1;
        elseif ~Temp_end1_A && ~Temp_end2_R
            E_AR_Cut=E_AR_Cut+1;
        elseif ~Temp_end1_R && ~Temp_end2_A
            E_AR_Cut=E_AR_Cut+1;
        elseif ~Temp_end1_R && ~Temp_end2_R
            E_R_Cut=E_R_Cut+1;
        elseif ~Temp_end1_M
            E_M1_Cut=E_M1_Cut+1;
        elseif ~Temp_end2_M
            E_M2_Cut=E_M2_Cut+1;
        else
            Ethers=Ethers+1;
        end
    end
    disp(q/Ptracts)
end

% convert tracts to mm (this is done in saving tracts function)
%for ti=1:t
%    Tracts{ti}(:,1)=(Tracts{ti}(:,1))*voxel_size(1);
%    Tracts{ti}(:,2)=(Tracts{ti}(:,2))*voxel_size(2);
%    Tracts{ti}(:,3)=(Tracts{ti}(:,3))*voxel_size(3);
%end


Tr='TRACK';
las='LAS';
INFO.id_string=Tr';
INFO.dim=Dimen(1:3)';
INFO.voxel_size=voxel_size';
origin(3)=0;
origin(2)=0;
origin(1)=0;
INFO.origin=origin';
INFO.n_scalar=0;
INFO.scalar_name=' ';
INFO.n_properties=0;
INFO.property_name=' ';
Vd=load_nii(V1_name,[],[],[],[],[],1);
INFO.vox_to_ras=Vd.hdr.hist.old_affine';
INFO.reserved='  ';
INFO.voxel_order=las';
INFO.pad2=las';
%INFO.img_orient_pat=[0; 0; 0; 0; 0; 0];
INFO.pad1=ls';
INFO.invert_x=0;
INFO.invert_y=0;
INFO.invert_z=0;
INFO.swap_xy=0;
INFO.swap_yz=0;
INFO.swap_zx=0;
INFO.n_count=t;
INFO.version=0;
INFO.hdr_size=1000;

% Write
timerun=toc;
Writetrk_rnh(INFO,Tracts,outputTRK)

disp(' ')
disp('Tractography general report ...')
disp('')
disp(['Initial Number of voxels Seeded = ',num2str(NSVoxels)])
disp(['Number of Seeds out of propagation ROI = ',num2str(NSVoxels-ANSVoxels)])
disp(['Number of potential tracts = ',num2str(Ptracts)])
disp('')
disp(['Num of tracts generated = ',num2str(t)])
disp(['Num of small tracts excluded = ',num2str(texcluded)])
disp('')
disp('Accepted tracts report ...')
disp(' ')
disp('Tracts with ends cutted for going out of propagation ROI:')
disp(T_R_Cut)
disp('Tracts with ends cutted for deviating to much:')
disp(T_A_Cut)
disp('Tracts with one end cutted for going out of ROI and the other for deviating:')
disp(T_AR_Cut)
disp('Tracts exceded maxima length in first way propagation')
disp(T_M1_Cut)
disp('Tracts exceded maxima length in second way propagation')
disp(T_M2_Cut)
disp('')
disp('Small excluded tracts report ...')
disp(' ')
disp('Tracts with ends cutted for going out of propagation ROI:')
disp(E_R_Cut)
disp('Tracts with ends cutted for deviating to much:')
disp(E_A_Cut)
disp('Tracts with one end cutted for going out of ROI and the other for deviating:')
disp(E_AR_Cut)
disp('Tracts exceded maxima length in first way propagation')
disp(E_M1_Cut)
disp('Tracts exceded maxima length in second way propagation')
disp(E_M2_Cut)
disp('')
disp('Time required for producing streamlines')
disp(timerun)

save([pth,sc,'Tractography', sc, outname, '_report_tdti'],...
'timerun',...
'NSVoxels',...
'ANSVoxels',...
'Ptracts',...
't',...
'texcluded',...
'T_A_Cut',...
'T_R_Cut',...
'T_AR_Cut',...
'T_M1_Cut',...
'T_M2_Cut',...
'E_A_Cut',...
'E_R_Cut',...
'E_AR_Cut',...
'E_M1_Cut',...
'E_M2_Cut',...
'Others',...
'Ethers');

save([pth,sc,'Tractography',sc, outname, 'M'], 'TractsM');