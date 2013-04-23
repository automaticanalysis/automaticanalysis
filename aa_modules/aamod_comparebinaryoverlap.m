function [aap,resp] = aamod_comparebinaryoverlap(aap, task, subjind)
% AAMOD_COMPAREBINARYOVERLAP Compare overlap of two binary images
%
% [aap,resp] = AAMOD_COMPAREBINARYOVERLAP(aap, task, subjind)
%
% input streams:   2 binary images
% output stream:   binaryoverlapmeasures
%
% Note that the first input stream is considered the "target" (if
% appropriate), the second the "test".
%
% Both are assumed to be in the same space and voxel-by-voxel matched.
%
% Overlap based largely on metrics from:
%
% Klein A, Andersson J, Ardekani BA, Ashburner J, Avants B, Chiang M-C,
% Christensen GE, Collins DL, Gee J, Hellier P, Song JH, Jenkinson M,
% Lepage C, Rueckert D, Thompson P, Vercauteren T, Woods RP, Mann JJ,
% Parsey RV (2009) Evaluation of 14 nonlinear deformation algorithms
% applied to human brain MRI registration. NeuroImage 46:786-802.


resp='';

switch task
    case 'domain'
        resp='subject';
    case 'description'
        resp='Compare overlap of two binary images';
    case 'doit'
                        
        targetStream = aap.tasklist.currenttask.inputstreams(1).stream{1};
        testStream = aap.tasklist.currenttask.inputstreams(1).stream{2};
        
        targetImg = aas_getfiles_bystream(aap, subjind, targetStream);
        testImg = aas_getfiles_bystream(aap, subjind, testStream);
        
        
        % read in images and reshape to be 1 by n vector
        Vtarget = spm_vol(targetImg);
        [Ytarget, XYZtarget] = spm_read_vols(Vtarget);
        Ytarget = reshape(Ytarget, 1, numel(Ytarget));
        
        Vtest = spm_vol(testImg);
        [Ytest, XYZtest] = spm_read_vols(Vtest);
        Ytest = reshape(Ytest, 1, numel(Ytest));
        
        
        % struecture for saving overlap measures
        overlapMeasures = [];
        
        
        % Target overlap
        overlapMeasures.targetOverlap = sum(Ytest.*Ytarget)./sum(Ytarget);
        
        % Mean overlap / Dice coefficient
        overlapMeasures.meanOverlap = 2 * (sum(Ytest.*Ytarget)./sum(Ytest + Ytarget));
        
        % Union overlap / Jaccard coefficient
        overlapMeasures.unionOverlap = (sum(Ytest.*Ytarget)./sum(Ytest + Ytarget));
        
        % False negative
        overlapMeasures.falseNegative = sum(Ytest < Ytarget); % 0 in test, 1 in target
        
        % False positive
        overlapMeasures.falsePostive = sum(Ytest > Ytarget); % 1 in test, 0 in target
        
        
        % get filepath and save overlapMeasures.mat                
        [pth, nm, ext] = fileparts(targetImg);
        outName = fullfile(pth, 'binaryOverlapMeasures.mat');
        save(outName, 'overlapMeasures');
        
        
       
        % describe output
        aap = aas_desc_outputs(aap, subjind, 'binaryoverlapmeasures', outName);        
end