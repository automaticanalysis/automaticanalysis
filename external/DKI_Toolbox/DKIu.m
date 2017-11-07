function DKIu
% United Diffusion Kurtosis Imaging toolbox
% This toolbox is for Diffusion Kurtosis processing - from essential
% pre-processing steps to the estimation of the kurtosis tensor and
% standard diffusion and kurtosis rotational invariant measures.
% This toolbox also include DKI biological modelling for estimates of
% axonal water fraction and also DKI based tractography
%
% Requirements: MATLAB and NIFTI_toolbox, images converted on NIFTI format,
% with bvalues and b-vectores extrated on two independent files (for the
% bvectors files, three lines have to correspond to three orientation 
% coordinates x,y and z
%
% References: Neto Henriques, R., Correia, M.M., Nunes, R.G., 
% Ferreira, H.A., 2015. Exploring the 3D Geometry of the Diffusion 
% Kurtosis Tensor - Impacts on the Development of Robust Tractography 
% Procedures and Novel Biomarkers. NeuroImage 111, 85-99. 
% doi.: 10.1016/j.neuroimage. 2015.02.004
%
% Developer: Rafael Neto Henriques
% Last update: 06/04/2015 

archstr = computer('arch');
st_comp=archstr(1:3);
if strcmp(st_comp,'win')
    sc='\';
else
    sc='/';
end

addpath(['..',sc,'NIFTI_toolbox'])

hfig=figure;

set(hfig,'color',[0.8 0.8 0.8],'name','DKIu',...
    'Units','normalized','NumberTitle','off','Toolbar',...
    'none','MenuBar','none','position',[0.3,0.2,0.2,0.6]);

uip=uipanel('FontSize',12,...
    'BackgroundColor',[0.93 0.93 0.93],'Units','Normalized',...
    'Position',[.025 .025 .95 .95],'HighlightColor','b',...
    'ForegroundColor','b');

% Logo
axes('position',[0.05  0.7672 0.9 0.1781]);
fill([0 0 8.2 8.2],[0 4 4 0],'black',[0 0 6],[0 4 0],'black',...
    [6 0 4 4 5 7 5 8],[0 4 4 2.7 4 4 2 0],[0 0 0.6], ...
    [8.2 8.2 8.9 8.9 9.3 9.3 10 10 9.6 8.6 8.2],...
    [0.3 2.2 2.2 1 1 2.2 2.2 0.3 0 0 0.3],[0 0 0.6])
xlim([0 10])
ylim([0 4])
axis off

% defaut strings
if strcmp(st_comp,'win')
    def1= 'C:\path\Subject1\DWI_brain.nii';
    def2='C:\path\Subject1\bvals.bval';
    def3='C:\path\Subject1\bvecs.bvec';
    def4='C:\path\Subject1\DWI_brain_mask.nii';
    def5='C:\path\Subject1\DWI_brain_GF150_olsDKI_DT.nii';
    def6='C:\path\Subject1\DWI_brain_GF150_olsDKI_KT.nii';
else
    def1='/path/Subject1/DWI_brain.nii';
    def2='/path/Subject1/bvals.bval';
    def3='/path/Subject1/bvecs.bvec';
    def4='/path/Subject1/DWI_brain_mask.nii';
    def5='/path/Subject1/DWI_brain_GF150_olsDKI_DT.nii';
    def6='/path/Subject1/DWI_brain_GF150_olsDKI_KT.nii';
end

% Button load data 
%uicontrol('Parent',uip,'Units','Normalized','Position',[0.05,0.70,0.9,0.06],...
%    'String','LOAD DATA','Callback',{@load_data});

% Button preprocessing Data
uicontrol('Parent',uip,'Units','Normalized','Position',[0.05,0.62,0.9,0.06],...
    'String','DWI Processing','Callback',{@DWI_Filter});

% Button DTI/DKI Fits
uicontrol('Parent',uip,'Units','Normalized','Position',[0.05, 0.54,0.9,0.06],...
    'String','DKI/DTI fit','Callback',{@fit_data});

% Button Standard DT metrics
uicontrol('Parent',uip,'Units','Normalized','Position',[0.05,0.46,0.9,0.06],...
    'String','Diffusion Standard Metrics','Callback',{@diff_stand_metrics});

% Button Standard KT metrics
uicontrol('Parent',uip,'Units','Normalized','Position',[0.05,0.38,0.9,0.06],...
    'String','Kurtosis Standard Metrics','Callback',{@kurt_stand_metrics});

% Button DKI biological modeling
uicontrol('Parent',uip,'Units','Normalized','Position',[0.05,0.30,0.9,0.06],...
    'String','DKI biological modeling','Callback',{@ad_metrics});

% Button Tractography
uicontrol('Parent',uip,'Units','Normalized','Position',[0.05,0.22,0.9,0.06],...
    'String','Tractography','Callback',{@dki_tract});

% % Butao Visualization
% uicontrol('Parent',uip,'Units','Normalized','Position',[0.05,0.1,0.9,0.06],...
%     'String','Visualization','Callback',{@Call_KurtVis});
% 
% Butao Simulation
uicontrol('Parent',uip,'Units','Normalized','Position',[0.05, 0.02,0.9,0.06],...
    'String','Simulation','Callback',{@Call_Rafael_to_simulate_your_problem});

%%
%     function load_data(source, eventdata)
%         p={'DWI file path:','bval file path:',...
%             'bvec file path:'};
%         def={def1;def2;def3};
%         r=inputdlg(p,'',1,def);
%         DWI_nii=r{1};
%         bval_bval=r{2};
%         bvec_bvec=r{3};
%         loaded_data=1;
%     end
%% Processing callback
    function DWI_Filter(source, eventdata)
        p={'DWI file path:'};
        def={def1};
        r=inputdlg(p,'',1,def);
        DWI_nii=r{1};
        Select_preprocessing_rnh(DWI_nii);
    end
%% DKI/DTI fitting callback
    function fit_data(source, eventdata)
        p={'DWI file path:','bval file path:','bvec file path:','Brain Mask path:'};
        def={def1;def2;def3;def4};
        r=inputdlg(p,'',1,def);
        DWI_nii=r{1};
        bval_bval=r{2};
        bvec_bvec=r{3};
        Mask_nii=r{4};
        Select_fitmethod_rnh(DWI_nii,Mask_nii,bval_bval,bvec_bvec);
    end
%% Diffusion measures callback
    function diff_stand_metrics(source, eventdata)
        p={'DT path:','Brain Mask path:'};
        def={def5;def4};
        r=inputdlg(p,'',1,def);
        DKI_DT=r{1};
        Mask_nii=r{2};
        fun_DTI_metrics(DKI_DT,Mask_nii);
    end
%% Kurtosis measures callback
    function kurt_stand_metrics(source, eventdata)
        p={'DT path:','KT path:','Brain Mask path:'};
        def={def5;def6;def4};
        r=inputdlg(p,'',1,def);
        DKI_DT=r{1};
        DKI_KT=r{2};
        Mask_nii=r{3};
        fun_DKI_metrics(DKI_DT,DKI_KT,Mask_nii);
    end
%% DKI Biological modeling callback
    function ad_metrics(source, eventdata)
        p={'DT path:','KT path:','Aligned fibers Mask:'};
        def={def5;def6;def4};
        r=inputdlg(p,'',1,def);
        DKI_DT=r{1};
        DKI_KT=r{2};
        Mask_svp_nii=r{3};
        fun_DKIadmetrics_rnh(DKI_DT,DKI_KT,Mask_svp_nii);
    end
%% DTI/DKI based tractography callback
    function dki_tract(source, eventdata)
        Select_tractography_rnh
    end
%%
    function Call_Rafael_to_simulate_your_problem(source, eventdata)
        DKIu_simulations
    end
%     function Call_KurtVis(source, eventdata)
%         KurtVis
%     end
end
