function [regionfile, labelfile, label, errorval]=getLabelMap(source)
%%% Written by Aaron Schultz (aschultz@martinos.org)
%%%
%%% Copyright (C) 2012,  Aaron P. Schultz
%%%
%%% Supported in part by the NIH funded Harvard Aging Brain Study (P01AG036694) and NIH R01-AG027435 
%%%
%%% This program is free software: you can redistribute it and/or modify
%%% it under the terms of the GNU General Public License as published by
%%% the Free Software Foundation, either version 3 of the License, or
%%% any later version.
%%% 
%%% This program is distributed in the hope that it will be useful,
%%% but WITHOUT ANY WARRANTY; without even the implied warranty of
%%% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%%% GNU General Public License for more details.

%% For proper functionallity with OrthoView, each possible source must be 
%   listed after the first if statement under regionfile. This is a unique 
%   case for when OrthoView is launched.
errorval=[];
label=[];
peak_nii_dir=fileparts(which('peak_nii.m'));
if nargin==0 || isempty(source)
    regionfile={'aal_MNI_V4' ...
        'Nitschke_Lab' ...
        'JHU_tracts' ...
        'JHU_whitematter' ...
        'Thalamus' ...
        'Talairach' ...
        'MNI' ...
        'HarvardOxford_cortex' ...
        'Juelich' ...
        'Cerebellum-flirt' ...
        'Cerebellum-fnirt' ...
        'Hammers' ...
        'Buckner7_loose' ...
        'Buckner7_tight' ...
        'Buckner17_loose' ...
        'Buckner17_tight' ...
        'LPBA40_flirt' ...
        'LPBA40_SPM5'};
else
    switch source
        case 'aal_MNI_V4'
            regionfile=[peak_nii_dir filesep 'aal_MNI_V4.img']; % From WFU_pickatlas
            labelfile=[peak_nii_dir filesep 'aal_MNI_V4_List.mat'];
            try 
                load(labelfile)
            catch 
                disp('Missing Atlas Files')
                errorval=['At least 1 Atlas File was not found: ' regionfile ' and'  labelfile];
                return
            end
        case 'Nitschke_Lab'
            regionfile=[peak_nii_dir filesep 'ControlabilityMask.nii']; % Provided by Deb Kerr
            labelfile=[peak_nii_dir filesep 'ControlabilityMask_List.mat'];
            try 
                load(labelfile)
            catch 
                disp('Missing Atlas Files')
                errorval=['At least 1 Atlas File was not found: ' regionfile ' and'  labelfile];
                return
            end
        case 'JHU_tracts'
            regionfile=[peak_nii_dir filesep 'JHU-ICBM-tracts-maxprob-thr0-1mm.nii']; % From FSL Atlas Files
            labelfile=[peak_nii_dir filesep 'JHU_tract_labels.mat'];
            try 
                load(labelfile)
            catch 
                disp('Missing Atlas Files')
                return
            end
        case 'JHU_whitematter'
            regionfile=[peak_nii_dir filesep 'JHU-WhiteMatter-labels-1mm.nii']; % From FSL Atlas Files
            labelfile=[peak_nii_dir filesep 'JHU_labels.mat'];
            try 
                load(labelfile)
            catch 
                disp('Missing Atlas Files')
                errorval=['At least 1 Atlas File was not found: ' regionfile ' and'  labelfile];
                return
            end
        case 'Thalamus'
            regionfile=[peak_nii_dir filesep 'Thalamus-maxprob-thr0-1mm.nii']; % From FSL Thalamus Atlas Files
            labelfile=[peak_nii_dir filesep 'Thalamus_labels.mat'];
            try 
                load(labelfile)
            catch 
                disp('Missing Atlas Files')
                errorval=['At least 1 Atlas File was not found: ' regionfile ' and'  labelfile];
                return
            end
        case 'Talairach'
            regionfile=[peak_nii_dir filesep 'Talairach-labels-1mm.nii']; % From FSL Talairach Atlas Files
            labelfile=[peak_nii_dir filesep 'Talairach_labels.mat'];
            try 
                load(labelfile)
            catch 
                disp('Missing Atlas Files')
                errorval=['At least 1 Atlas File was not found: ' regionfile ' and'  labelfile];
                return
            end
        case 'MNI'
            regionfile=[peak_nii_dir filesep 'MNI-maxprob-thr0-1mm.nii']; % From FSL MNI Atlas Files
            labelfile=[peak_nii_dir filesep 'MNI_labels.mat'];
            try 
                load(labelfile)
            catch 
                disp('Missing Atlas Files')
                errorval=['At least 1 Atlas File was not found: ' regionfile ' and'  labelfile];
                return
            end
        case 'HarvardOxford_cortex'
            regionfile=[peak_nii_dir filesep 'HarvardOxford-cort-maxprob-thr0-1mm.nii']; % From FSL Atlas Files
            labelfile=[peak_nii_dir filesep 'HarvardOxford_cortical_labels'];
            try 
                load(labelfile)
            catch 
                disp('Missing Atlas Files')
                errorval=['At least 1 Atlas File was not found: ' regionfile ' and'  labelfile];
                return
            end
            %case 'HarvardOxford_subcortical' % Labels do not match image
            %    regionfile=[peak_nii_dir filesep 'HarvardOxford-sub-maxprob-thr0-1mm.nii']; % From FSL Atlas Files
            %    labelfile=[peak_nii_dir filesep 'HarvardOxford_subcortical_labels']);
            %    load(labelfile)
        case 'Juelich'
            regionfile=[peak_nii_dir filesep 'Juelich-maxprob-thr0-1mm.nii']; % From FSL Atlas Files
            labelfile=[peak_nii_dir filesep 'Juelich_labels.mat'];
            try 
                load(labelfile)
            catch 
                disp('Missing Atlas Files')
                errorval=['At least 1 Atlas File was not found: ' regionfile ' and'  labelfile];
                return
            end
        case 'Cerebellum-flirt'
            regionfile=[peak_nii_dir filesep 'Cerebellum-MNIflirt-maxprob-thr0-1mm.nii']; % From FSL Atlas Files
            labelfile=[peak_nii_dir filesep 'Cerebellum_labels.mat'];
            try 
                load(labelfile)
            catch 
                disp('Missing Atlas Files')
                errorval=['At least 1 Atlas File was not found: ' regionfile ' and'  labelfile];
                return
            end
        case 'Cerebellum-fnirt'
            regionfile=[peak_nii_dir filesep 'Cerebellum-MNIfnirt-maxprob-thr0-1mm.nii']; % From FSL Atlas Files
            labelfile=[peak_nii_dir filesep 'Cerebellum_labels.mat'];
            try 
                load(labelfile)
            catch 
                disp('Missing Atlas Files')
                errorval=['At least 1 Atlas File was not found: ' regionfile ' and'  labelfile];
                return
            end
        case 'Hammers'
            if exist([fileparts(which('spm')) filesep 'toolbox' filesep 'HammersAtlas' filesep 'Hammers_mith_atlas_n30r83_SPM5.img'],'file') && exist([fileparts(which('spm')) filesep 'toolbox' filesep 'HammersAtlas' filesep 'Hammers_mith_atlas_n30r83_labels.mat'],'file')
                regionfile=[fileparts(which('spm')) filesep 'toolbox' filesep 'HammersAtlas' filesep 'Hammers_mith_atlas_n30r83_SPM5.img'];
                labelfile=[peak_nii_dir filesep 'Hammers_mith_atlas_n30r83_labels.mat'];
                try 
                load(labelfile)
            catch 
                disp('Missing Atlas Files')
                errorval=['At least 1 Atlas File was not found: ' regionfile ' and'  labelfile];
                return
            end
            else
                errorval='The HammersAtlas is not available. Please contact Alex Hammers (alexander.hammers@fondation-neurodis.org) for the atlas';
                disp('===========')
                disp('IMPORTANT: Once you have downloaded the Atlas files, please put them into a directory called HammersAtlas in the toolbox directory of SPM');
                disp('===========')
                label.source=source;
                return
            end
        case 'Buckner7_loose'
            regionfile=[peak_nii_dir filesep 'Buckner7_loose.nii']; % Derived from http://surfer.nmr.mgh.harvard.edu/fswiki/CerebellumParcellation_Buckner2011 & http://surfer.nmr.mgh.harvard.edu/fswiki/CorticalParcellation_Yeo2011
            labelfile=[peak_nii_dir filesep 'BucknerLabels.mat'];
            try 
                load(labelfile)
            catch 
                disp('Missing Atlas Files')
                errorval=['At least 1 Atlas File was not found: ' regionfile ' and'  labelfile];
                return
            end
        case 'Buckner7_tight'
            regionfile=[peak_nii_dir filesep 'Buckner7_tight.nii']; % Derived from http://surfer.nmr.mgh.harvard.edu/fswiki/CerebellumParcellation_Buckner2011 & http://surfer.nmr.mgh.harvard.edu/fswiki/CorticalParcellation_Yeo2011
            labelfile=[peak_nii_dir filesep 'BucknerLabels.mat'];
            try 
                load(labelfile)
            catch 
                disp('Missing Atlas Files')
                errorval=['At least 1 Atlas File was not found: ' regionfile ' and'  labelfile];
                return
            end
        case 'Buckner17_loose'
            regionfile=[peak_nii_dir filesep 'Buckner17_loose.nii']; % Derived from http://surfer.nmr.mgh.harvard.edu/fswiki/CerebellumParcellation_Buckner2011 & http://surfer.nmr.mgh.harvard.edu/fswiki/CorticalParcellation_Yeo2011
            labelfile=[peak_nii_dir filesep 'BucknerLabels.mat'];
            try 
                load(labelfile)
            catch 
                disp('Missing Atlas Files')
                errorval=['At least 1 Atlas File was not found: ' regionfile ' and'  labelfile];
                return
            end
        case 'Buckner17_tight'
            regionfile=[peak_nii_dir filesep 'Buckner7_tight.nii']; % Derived from http://surfer.nmr.mgh.harvard.edu/fswiki/CerebellumParcellation_Buckner2011 & http://surfer.nmr.mgh.harvard.edu/fswiki/CorticalParcellation_Yeo2011
            labelfile=[peak_nii_dir filesep 'BucknerLabels.mat'];
            try 
                load(labelfile)
            catch 
                disp('Missing Atlas Files')
                errorval=['At least 1 Atlas File was not found: ' regionfile ' and'  labelfile];
                return
            end
            %case 'Custom'
            %    regionfile=['ROI image file'];
            %    labelfile=['ROIlist mat-file']);
        case 'LPBA40_flirt'
            if exist([fileparts(which('spm')) filesep 'toolbox' filesep 'LBPA40_flirt' filesep 'maxprob' filesep 'lpba40.flirt.avg152T1_brain.label.nii'],'file') && exist([fileparts(which('spm')) filesep 'toolbox' filesep 'HammersAtlas' filesep 'Hammers_mith_atlas_n30r83_labels.mat'],'file')
                regionfile=[fileparts(which('spm')) filesep 'toolbox' filesep 'LBPA40_flirt' filesep 'maxprob' filesep 'lpba40.flirt.avg152T1_brain.label.nii'];
                labelfile=[peak_nii_dir filesep 'lbpa40_labels.mat'];
                try 
                load(labelfile)
            catch 
                disp('Missing Atlas Files')
                errorval=['At least 1 Atlas File was not found: ' regionfile ' and'  labelfile];
                return
            end
            elseif exist([fileparts(which('spm')) filesep 'toolbox' filesep 'LBPA40_flirt' filesep 'maxprob' filesep 'lpba40.flirt.avg152T1_brain.label.nii.gz'],'file') && exist([fileparts(which('spm')) filesep 'toolbox' filesep 'HammersAtlas' filesep 'Hammers_mith_atlas_n30r83_labels.mat'],'file')
                gunzip([fileparts(which('spm')) filesep 'toolbox' filesep 'LBPA40_flirt' filesep 'maxprob' filesep 'lpba40.flirt.avg152T1_brain.label.nii.gz'])
                regionfile=[fileparts(which('spm')) filesep 'toolbox' filesep 'LBPA40_flirt' filesep 'maxprob' filesep 'lpba40.flirt.avg152T1_brain.label.nii'];
                labelfile=[peak_nii_dir filesep 'lbpa40_labels.mat'];
                try 
                load(labelfile)
            catch 
                disp('Missing Atlas Files')
                errorval=['At least 1 Atlas File was not found: ' regionfile ' and'  labelfile];
                return
            end
            else
                errorval='The LBPA40 is not available. Please download the atlas from http://www.loni.ucla.edu/Protocols/LPBA40';
                disp('===========')
                disp('IMPORTANT: Once you have downloaded the Atlas files, please put them into a directory called LBPA40_flirt in the toolbox directory of SPM');
                disp('===========')
                label.source=source;
                return
            end
        case 'LPBA40_SPM5'
            if exist([fileparts(which('spm')) filesep 'toolbox' filesep 'LBPA40_SPM5' filesep 'maxprob' filesep 'lpba40.SPM5.avg152T1.label.nii'],'file') && exist([fileparts(which('spm')) filesep 'toolbox' filesep 'HammersAtlas' filesep 'Hammers_mith_atlas_n30r83_labels.mat'],'file')
                regionfile=[fileparts(which('spm')) filesep 'toolbox' filesep 'LBPA40_SPM5' filesep 'maxprob' filesep 'lpba40.SPM5.avg152T1.label.nii'];
                labelfile=[peak_nii_dir filesep 'lbpa40_labels.mat'];
                try 
                load(labelfile)
            catch 
                disp('Missing Atlas Files')
                errorval=['At least 1 Atlas File was not found: ' regionfile ' and'  labelfile];
                return
            end
            elseif exist([fileparts(which('spm')) filesep 'toolbox' filesep 'LBPA40_SPM5' filesep 'maxprob' filesep 'lpba40.SPM5.avg152T1.label.nii.gz'],'file') && exist([fileparts(which('spm')) filesep 'toolbox' filesep 'HammersAtlas' filesep 'Hammers_mith_atlas_n30r83_labels.mat'],'file')
                gunzip([fileparts(which('spm')) filesep 'toolbox' filesep 'LBPA40_SPM5' filesep 'maxprob' filesep 'lpba40.SPM5.avg152T1.label.nii.gz'])
                regionfile=[fileparts(which('spm')) filesep 'toolbox' filesep 'LBPA40_SPM5' filesep 'maxprob' filesep 'lpba40.SPM5.avg152T1.label.nii'];
                labelfile=[peak_nii_dir filesep 'lbpa40_labels.mat'];
                try 
                load(labelfile)
            catch 
                disp('Missing Atlas Files')
                errorval=['At least 1 Atlas File was not found: ' regionfile ' and'  labelfile];
                return
            end
            else
                errorval='The LBPA40 is not available. Please download the atlas from http://www.loni.ucla.edu/Protocols/LPBA40';
                disp('===========')
                disp('IMPORTANT: Once you have downloaded the Atlas files, please put them into a directory called LBPA40_SPM5 in the toolbox directory of SPM');
                disp('===========')
                label.source=source;
                return
            end
        otherwise
            regionfile=[];
            labelfile=[];
            errorval='Source name did not match any atlases';
            return            
    end
    
    label.source=source;
    label.rf.hdr=spm_vol(regionfile);
    [label.rf.img,label.rf.XYZmm]=spm_read_vols(label.rf.hdr);
    label.ROI=ROI;
    label.ROInames={ROI.Nom_C}; %ROInames is taken from ROI, where ROI is structure
    if ~isfield(ROI,'ID')
        display('ERROR: ROI must be a structure with an ID field that is the ID of values in region file')
        return
    end
end