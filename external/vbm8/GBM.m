function [src,cls,mask,TSEG,SLAB,opt] = GBM(src,cls,res,opt)
% SKULL-STRIPPING USING GRAPH-CUT
% _________________________________________________________________________
%
% scr (single)    bias corrected T1 image after noise reduction
% cls (cell)      6 uint8 tissue class images (GM,WM,CSF,.,.,head)
% res (struct)    
% opt (struct)    main parameter
%  .RSS           remove sinus sagittalis 
%                 controls the detection (has to be greater or equal one) 
%                 (RSS=1 more tissue, RSS=2.0 medium tissue (default), RSS>2 less tissue)  
%  .BVH           remove high intensity blood vessels
%  .BVN           remove low  intensity blood vessels
%  .FILL          detect ventricle and deep GM structures
%                 create segment map with filled ventrile and deep GM structures
%
%
% controll parameter:
% RSS  = remove sinus sagittalis 
%        controls the detection (has to be greater or equal one)
%        (RSS=1 more tissue, RSS=1.5 medium tissue (default), RSS>=2 less tissue)
%
% RBV  = remove blood vessels
% FILL = create segment map with filled ventrile and deep GM structures
% _________________________________________________________________________
% Robert Dahnke 2011_01
% Center of Neuroimaging 
% University Jena
% $Id: GBM.m 423 2011-07-07 12:15:48Z gaser $

% CHECK INPUT
% _________________________________________________________________________

  if ~exist('opt','var'), opt=struct(); end

  def.RSS    = 2.0; % RSS controls the detection has to be greater or equal one (RSS=1 more tissue, RSS=1.5 medium tissue (default), RSS>=2 less tissue)
  def.BVH    = 1;   % remove high intensity blood vessels
  def.BVN    = 0;   % remove low  intensity blood vessels
  if nargout > 3     % detection only if TSEG and SLAB are defined as output
      def.FILL   = 1;   % detect ventricle and deep GM structures
  else
      def.FILL   = 0;   % don't detect ventricle and deep GM structures
  end
  def.verb   = 1;   % display process 
  
  opt = checkinopt(opt,def);
  
  vx_vol      = sqrt(sum(res.image(1).mat(1:3,1:3).^2));
  scale_morph = 1/mean(vx_vol);
  
% first simple correction of classes (not ready!)
% _________________________________________________________________________
  stime = clock;
  tmp = sum(single(cat(4,cls{1},cls{2},cls{3},cls{4},cls{5},cls{6})),4);
  tmp(tmp>255) = 255;
  cls{6}=cls{6} + uint8(255-tmp); % if the sum of the classe is not 100% add the rest to the background 
  for i=1:6, cls{i}=uint8(median3(single(cls{i}),cls{i}>0,true(size(src)))); end              % noise reduction  

  % 2) if the head is to thick, SPM detect GM outside the brain and we have to remove it (only one GM Segment)
  sum_cls12 = single(cls{1}) + single(cls{2});
  M = single(cg_morph_vol(cg_morph_vol(sum_cls12>64,'open',3),'labopen')); M(sum_cls12<64)=-inf;   
  M = cg_morph_vol(cg_morph_vol(uint8(down_cut(M,sum_cls12/255*3,0.1,vx_vol)==1),'labopen'),'dilate'); 
  M = cg_morph_vol(cg_morph_vol(cg_morph_vol(uint8((single(cls{1}.*M)+single(cls{2}))>64),'open',2),'labopen'),'dilate',3); cls{1}=cls{1} .* M; 
  % 1) there should be only one WM Segment
  M = cg_morph_vol(cg_morph_vol(uint8(cls{2}>64),'labopen'),'dilate'); cls{2}=cls{2} .* M; %cls{1}=cls{1} + cls{2} .* (1-M); % 1.
  % 2) big ventricle can lead to background problems, i.e. IXI175
  % 3) big ventricle can lead to problems with the probability map - resulting in WM within the ventricle, i.e. 4397
  
  
    
% Inital atlas map SLAB and create the TSEG map a tissue-based scaled 
% Version of the src image (BG=0, CSF=1, GM=2, WM=3, Other=4). 
% _________________________________________________________________________
  clstp=[0 1 0 11 11 11]; clsth=[1 0.9 1 0.9 0.6 0.2]; 
  sum_cls12 = single(cls{1}) + single(cls{2});
  SLAB = zeros(size(src),'single'); for c=1:6, SLAB(cls{c}>255*clsth(c))=clstp(c); end; 
  mask = cg_morph_vol(cg_morph_vol(cg_morph_vol(sum_cls12/255*3>0.5,'open',round(2*scale_morph))==1,'dilate',round(2*scale_morph)),'labopen');
  SLAB((cls{4}>255*0.8 | cls{5}>255*0.5 | cls{6}>255*0.5) & mask==0)=11; clear mask;

  % inital speed map for graph-cut (BG=0, CSF=1, GM=2, WM=3, Other=4) 
  % because we adjust the intensities noise filter can used again 
  % TODO: add sub-classes and classes for other tissues (BV=4, HD=5)
  TSEG = create_TSEG(src,cls,vx_vol); 
  %sanlmMex(TSEG,3,1);
    
  
% REMOVE SINUS SAGITALIS
% _________________________________________________________________________
% Remove sinus sagittalis by analysis of the range from 1.25 to 1.75
% between GM and CSF.
  cls4=0;
  if opt.RSS>0
    sum_cls123 = sum_cls12 + single(cls{3});
    sum_cls456 = single(cls{4}) + single(cls{5}) + single(cls{6});
    M = single(sum_cls123>sum_cls456); 
    clear sum_cls456
    M( cg_morph_vol( ((single(cls{1})/3)>cls{3} & cls{3}<64) | cls{2}>0,'open',1)==1 & M==1)=-1; D1=eikonal3(M,vx_vol);
    M = single(cg_morph_vol(cg_morph_vol(D1/opt.RSS>TSEG,'close'),'erode')); 
    M(M==0 & cg_morph_vol(cg_morph_vol( (sum_cls123/255<0.25) | (sum_cls12/255>0.75) ,'open'),'close')==1)=-inf; 
    clear sum_cls123
    M = cg_morph_vol(down_cut(M,TSEG,0.1,vx_vol),'close'); 

    % removing of inbrain structures (simple closing will also remove small structure that we are interested in 
    D2 = single(cg_morph_vol(M==1 | SLAB==11,'labclose')); D2(cg_morph_vol(cg_morph_vol(D2,'open',6*scale_morph),'labopen')==1)=-1; 
    D2 = eikonal3(D2,vx_vol); D2=median3(D2,D2>0); 
    S  = false(size(TSEG)); S((cls{4}<255*0.25 & cls{5}<255*0.25 & cls{6}<255*0.1) & M==1 & D2<10 & D1>TSEG)=1; clear D1 M;
  else
    D2 = zeros(size(SLAB),'single');
    S  = zeros(size(SLAB),'single');
  end
  
  
    
% HIGH INTENSITY BLOOD VESSELS 
% _________________________________________________________________________
  BV=(cg_morph_vol(cg_morph_vol(SLAB==11,'open',3),'labopen')==0 & SLAB==11) | (TSEG>3.25 & SLAB<=1); 

  
% MAJOR ALIGNMENT       
% _________________________________________________________________________
% align other voxels and smooth the map (median filter)
  SLAB((SLAB==0 & TSEG<1.75) | BV)=-inf; SLAB = down_cut(SLAB,TSEG,0.0,vx_vol); SLAB(SLAB==-inf)=0; clear TSEGC;  % alignment without sinus sagittalis
  M=cg_morph_vol(D2<=inf,'labopen'); SLAB(S>0 & M==1 & SLAB~=1)=11; SLAB(SLAB==11 & M==0)=7;                 % add sinus sagittalis anc other BV
  if opt.BVH, SLAB(BV)=7; else SLAB(BV)=0; end
  SLAB = down_cut(SLAB,TSEG,0.01,vx_vol);                                                        % align only neighbors with lower tissue intentsity
  SLAB = down_cut(SLAB,TSEG,4.00,vx_vol);                                                        % align all remaining voxels
  SLAB = median3(SLAB,TSEG>0.5);                                                                   % smoothing
  clear S; opt.time.GBM8GBM = dp('GBM',opt,stime); stime = clock;



% NORMAL INTENSITY BLOOD VESSELS
% _________________________________________________________________________
% Calculate GM-WM tissue thickness to remove elements that are thinner than t mm.  
% TODO: - use D2 to descide between VB and HD ... 
%       - if you find BVs and want to use them for further analaysis
%         (to reduce fMRI errors), you need to downcut them
  if opt.BVN
    T = single( (sum_cls12/255<0.1) | SLAB==11);
    D = vbdist(T,true(size(TSEG)),vx_vol);
    T = thickness_V0x59_TIM(T+2,D,zeros(size(D),'single')); 
    T = median3(T,TSEG>0.5); 
    t = 2; M = (T<=t) & (SLAB==1) & ((sum_cls12+single(cls{4})+single(cls{5}))/255>0.5); SLAB(cg_morph_vol(M | SLAB==11,'erode',1)==1 & SLAB~=11)=0;
    SLAB = down_cut(SLAB,TSEG,-0.1,vx_vol);
    SLAB(median3(single(cg_morph_vol(cg_morph_vol(SLAB==1 & T>0 & T<=t/2 & sum_cls12/255>0.25,'dilate'),'erode',1)))==1)=7;
    TSEGC=TSEG; TSEGC(sum_cls12/255<0.1 | T>2*t)=-inf; SLAB(cg_morph_vol(down_cut(single(SLAB==7),TSEGC,0.01,vx_vol),'close')==1 & SLAB==1)=7;
    clear D T;
    opt.time.GBM8BV = dp('|BV',opt,stime); stime = clock;  
  end

  clear sum_cls12
    
% BRAINMASK
% _________________________________________________________________________
% Use the detected region (SLAB==1), remove holes and islands and add some
% space around the the major segment.
  dc = 2;
  mask  = (SLAB==1); % | SLAB==7); % | (D2>5);                                 % initial mask, D2 is only available if RSS>0)
  mask  = cg_morph_vol(mask,'close',1);
  mask  = cg_morph_vol(mask,'labclose');                                   % closing holes
  mask  = cg_morph_vol(mask,'labopen');                                    % remove islands
  % add some tissue arround the mask, but only if the intensity is lower and if there is no BV
  D = vbdist(single(mask==1),true(size(SLAB)),vx_vol);
  C = vbdist(single(D>dc),true(size(SLAB)),vx_vol);
  SLAB(D>0 & D<=3.0 & SLAB~=7)=0;
  SLAB  = down_cut(SLAB,TSEG,0.0,vx_vol);                                % align only neighbors with lower tissue intentsity
  SLAB  = down_cut(SLAB,TSEG,4.0,vx_vol);                                % align all remaining voxels
  SLAB2 = median3(SLAB); SLAB(SLAB~=7)=SLAB2(SLAB~=7); clear SLAB2;        % don't filter blood vessels!
  mask  = median3(single(mask | SLAB==1 | (D<=1.5) | C>(dc*1.25)))>0.5;    % final brain mask
  clear ROI D clsth clstp dth;
  opt.time.GBM8M = dp('|M',opt,stime); stime = clock;
  
  

  if opt.FILL
  % VENTRICLE DETECTION
  % _____________________________________________________________________
  % Detect ventricle by creating sulcal depth. Voxel that are higher than
  % 90% and voxel that are lower than 10% of the maximum are used as seed
  % regions. Growing is regularised by images intensity to avoid
  % overflowing of the small boundaries. Another problems ist that, we
  % have to find a good threshold to calculate the eikonal distance. If
  % the threshold is too low, we may not able to get to the ventricle.
  % If the threshold is too high we remove important boundaries and we get
  % from the wrong side (normaly next to the hindbrain) access to the
  % ventricle and the values are too low, leading to a lot of errors.
  % Small ventricle can lead to detection problmes because the are too
  % thin to allow greater depth.
  % TODO: - something to check the result (relation between ventricle and non-ventricle? - no?, to low max depth?)
  %       = growing process looks good, but the detection is horrible
  %       - detection: - labeling of csf-segment 
  %       = for the distance thresholds you has to use the maximum distnaces
    foundV=0; TSEGth=1.5:0.5:2.5; i=0; Dmax=0;
    DM=vbdist(single(cg_morph_vol(cg_morph_vol(SLAB==11 | SLAB==8,'open',2),'labopen')),true(size(SLAB)),vx_vol); DMmax=max(DM(DM<inf));
    while foundV==0 && i<numel(TSEGth) && (Dmax<1.5*DMmax);
      i=i+1; substep=0:1; j=0; allCSF=0;
      while allCSF==0 && j<numel(substep)
        j=j+1; TSEGa=substep(j)*(TSEGth(2)-TSEGth(1))/numel(substep);
        M=single(TSEG<TSEGth(i)+TSEGa); M(DM==0)=-1; sumCSF=sum(M(:)==1); 
        D=eikonal3(M,vx_vol); Dmax=max(D(D<inf)); sumINF=sum(D(M(:)==1)==inf);                            % calculate sulcal depth
        allCSF = (sumINF/sumCSF)<0.01;
      end
      M=single(D<inf & (D>(0.9*Dmax) | (D>(0.8*1.5*DMmax))) & TSEG<(TSEGth(i)+TSEGa));                  % set ventricular starting region
      foundV=sum(M(:)==1); 
    end
    if foundV>0
      % find the two greatest objects - you have do do more than that...
      % 1) there are maybe less elements... 
      % 2) may both ventricle are one...

      M(cg_morph_vol(TSEG>=1,'close')==1 & M~=1)=-inf; M=down_cut(M,max(0,3-TSEG),0.25,vx_vol);        % growing only for ventricular region 
      M(cg_morph_vol(SLAB==11 | SLAB==8 | (D>0 & D<(0.1*Dmax)),'labopen')==1)=11;                        % set non-ventricular starting region
      for ti=1:0.25:2, M(M==-inf & TSEG<=ti)=0; M=down_cut(M,max(0,2-TSEG),0.25,vx_vol); end             % levelwise growing of both regions
      M1=cg_morph_vol(cg_morph_vol(M==1,'open'),'labopen'); M(M1==0 & M==1 )=0; clear M1;                 % only the biggest object
      SLAB(M==1)=15; SLAB(M==11 & SLAB==1)=21;                                                            % add ventricle and non-ventricle to the atlas
      opt.time.GBM8V = dp('|V',opt,stime); stime = clock;

      

    % BASAL STRUCTURES
    % _____________________________________________________________________
    % Detect inner GM structures using the distance to ventricles DV, TSEG 
    % intensity (class between GM and WM and region growing.
    % Mainly controlled by dbg: 
    % (1) max distance to the ventricle for the region for the growing process
    % (2) max distance to the ventricle for the start region
    % (3) min distance to the non-ventricle region for the start region an region growing
    % (4) first  max distance to the start region for the growing process 
    % (5) second max distance to the start region for the growing process
    % TODO: - strong relation to the ventricle
    %       = division into sup-structures
      dbg = [30 15 8 6 10]; 
      DV  = vbdist(single(SLAB==15),true(size(SLAB)),vx_vol); 
      DNV = vbdist(single(SLAB==21),true(size(SLAB)),vx_vol);
      M=-inf(size(SLAB),'single'); 
      M(cg_morph_vol(cg_morph_vol((SLAB==1 & DV<dbg(1) & SLAB~=21 & TSEG>2.0 & TSEG<2.75) | SLAB==15,'close'),'open',1)==1)=0; % intensity region
      M(M==0 & cg_morph_vol(SLAB~=15,'erode',2)==1 & DV<dbg(2) & DNV>dbg(3))=1;                                                % starting region 
      D=vbdist(single(M>0),true(size(SLAB)),vx_vol); % dist from start region (may replace by eikonal distance)
      % region growing (but only near the starting region)
      M(D>dbg(4) | SLAB==15 | TSEG<2.25)=-inf;                            M=down_cut(M,        TSEG ,0,vx_vol); % growing to the GM
      M(D<dbg(4) & SLAB==1  & TSEG<2.75 & M==-inf)=0;                     M=down_cut(M,max(0,3-TSEG),0,vx_vol); % growing to the WM
      M(D<dbg(5) & SLAB==1  & TSEG<2.75 & TSEG>1.00 & M==-inf & DV<10)=0; M=down_cut(M,        TSEG ,0,vx_vol); % growing to the ventricle
      M=median3(M); M(cg_morph_vol(M==-inf,'erode',2)==0 & M==-inf)=0; M=median3(M);                              % smoothing
      SLAB(SLAB==1 & M==1 & TSEG<=2.75)=5;                                                                        % add ventricle to the atlas 
      opt.time.GBM8BG = dp('|BG',opt,stime); stime = clock;


    % CEREBELLUM
    % _____________________________________________________________________
    % ... find it and than you can decide beween left and right with the BG
    % structures


    % HEMISPHERES
    % _____________________________________________________________________
    % TODO: - detection: - use Basalganlia and the center of mass to find the hemispheres
    %                    - to descide which one is left/right use the mat data 
    %                    - 
    else
       if opt.verb>0, fprintf('|VD: failed     '); end
    end
  end

function T=create_TSEG(src,cls,vx_vol)
% _________________________________________________________________________
% Scales the intensity of the T1 images (T1f) similar to the segment image
% SEGf. Furthermore, only segment data will used for tisses above 2 (GM).
% if no T1 file is given it try to find the modulated T1 file of the VBM
% segmention.
% If there is no T1 image a waring will displayed and TSEG is given by SEG.
% 
%   T=create_TSEG(src,cls,opt)
% _________________________________________________________________________
% Center of Neuroimaging, University Jena, Germany
% Robert Dahnke
% 2011/01

  % Scale the T1 data based on the SEG information (tissue mean values)
  th=[0.5 0.5 0.5]; minsrc=min(src(:)); maxsrc=max(src(:)); dth=[1 4 20 0 0 0];
  srcs=smooth3(src); % smooth the src to get more stable histograms
  % we can't trust the CSF region outside the brainm that why we only want to take values from voxel with a creater distance to the head
  sum_cls123 = single(cls{1}) + single(cls{2}) + single(cls{3});
  sum_cls456 = single(cls{4}) + single(cls{5}) + single(cls{6});
  D = vbdist(single(cg_morph_vol(sum_cls456>sum_cls123,'close',2)),true(size(src)),vx_vol);
  clear sum_cls123 sum_cls456
  
  [hst,hsti]=hist(srcs(cls{1}(:)/255>th(1) & D(:)>dth(1)),linspace(minsrc,maxsrc,1000));               [tmp,maxi]=max(hst); t(2)=hsti(maxi); % GM
  [hst,hsti]=hist(srcs(cls{2}(:)/255>th(2) & D(:)>dth(2) & src(:)>t(2)),linspace(minsrc,maxsrc,1000)); [tmp,maxi]=max(hst); t(1)=hsti(maxi); % WM
  [hst,hsti]=hist(srcs(cls{3}(:)/255>th(3) & D(:)>dth(3) & src(:)<t(2)),linspace(minsrc,maxsrc,1000)); [tmp,maxi]=max(hst); t(3)=hsti(maxi); % CSF
  clear tmp hst hsti maxi;
  
  T = single(src);
  M = T<=t(3);        T(M) =   (T(M) - min(T(:))) / max(eps,(t(3)-min(T(:))));    % correct CSF
  M = T>1 & T<=t(2);  T(M) = 1+(T(M) - t(3))      / (t(2)-t(3));                  % correct GM
  M = T>2 & T<=t(1);  T(M) = 2+(T(M) - t(2))      / (t(1)-t(2));                  % correct WM
  M = T>3;            T(M) = 3+(T(M) - t(1))      / (max(T(:))-t(1));             % correct WM+
  
  % If someone has thresholded the image and there is no difference between
  % background and CSF. The other threshold (WM) should be no problem.
  if t(3)==minsrc
    M=T<1 & (cls{3}>0 | cls{2}>0); T(M) = max((single(cls{3}(M)) + single(cls{1}(M)))/255*0.75,T(M));          % avoid background in the brain
    M=T>0 & (single(cls{6})>(single(cls{1})+single(cls{2})+single(cls{3})+single(cls{4})+single(cls{5}))); T(M) = min(1-single(cls{6}(M)/255),T(M)); % clear background
    M=smooth3(T);T(T<1)=M(T<1);
  end
  
return