function [aap,resp] = aamod_comparecontinuousoverlap(aap, task, subjind)
% AAMOD_COMPARECONTINUOUSOVERLAP Compare overlap of two binary images
%
% [aap,resp] = AAMOD_COMPARECONTINUOUSOVERLAP(aap, task, subjind)
%
% input streams:   1 binary image, 1 continuous
% output stream:   continuousoverlapmeasures
%
% Note that the first input stream is considered the "target" and should be
% binary, the second the "test" (and continuous).
%
% Both are assumed to be in the same space and voxel-by-voxel matched.
%
% Overlap based largely on metrics from:
% 
% Crum WR, Camara O, Rueckert D, Bhatia KK, Jenkinson M, Hill DLG (2005)
% Generalised overlap measures for assessment of pairwise and groupwise
% image registration and segmentation. MICCAI 3749:99-106.
%
% Klein A, Andersson J, Ardekani BA, Ashburner J, Avants B, Chiang M-C,
% Christensen GE, Collins DL, Gee J, Hellier P, Song JH, Jenkinson M,
% Lepage C, Rueckert D, Thompson P, Vercauteren T, Woods RP, Mann JJ,
% Parsey RV (2009) Evaluation of 14 nonlinear deformation algorithms
% applied to human brain MRI registration. NeuroImage 46:786-802.
% 


resp='';

switch task
    case 'domain'
        resp='subject';
    case 'description'
        resp='Compare overlap of continuous image with a binary image';
    case 'doit'
                        
        targetStream = aap.tasklist.currenttask.inputstreams(1).stream{1};
        testStream = aap.tasklist.currenttask.inputstreams(1).stream{2};        
        
        targetImg = aas_getfiles_bystream(aap, subjind, targetStream);
        testImg = aas_getfiles_bystream(aap, subjind, testStream);
        
        
        % read in images and reshape to be 1 by n vector
        Vtarget = spm_vol(targetImg);
        [Ytarget, XYZtarget] = spm_read_vols(Vtarget);
        Ytarget = reshape(Ytarget, 1, numel(Ytarget));
        targetTrueInd = find(Ytarget==1);
        targetFalseInd = find(Ytarget==0); 
        
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
        
        imageDiff = Ytest-Ytarget;
        
        % MSE
        overlapMeasures.MSE = mean(imageDiff.^2);
        
        % False negative (where imageDiff is negative)
        whereNegative = find(imageDiff<0);
        overlapMeasures.falseNegative = sum(Ytest(whereNegative)); % 0 in test, 1 in target
        
        % False positive
        wherePositive = find(imageDiff>0);
        overlapMeasures.falsePostive = sum(Ytest(wherePositive)); % 1 in test, 0 in target
        
        % d'
        % d' = z(H) - z(F)
        hits = (sum(Ytest(targetTrueInd))/sum(targetTrueInd)) / length(targetTrueInd);
        falseAlarms = (sum(Ytest(targetFalseInd))/sum(targetFalseInd)) / length(targetFalseInd);
        
        if hits==1
            hits = length(Ytarget)/(length(Ytarget)+1);
        end
        
        if falseAlarms==0
            falseAlarms = 1/(length(Ytarget)+1);
        end
        
        overlapMeasures.dprime = norminv(hits, 0, 1) - norminv(falseAlarms, 0, 1);
        
        
        
        % get filepath and save overlapMeasures.mat                
        [pth, nm, ext] = fileparts(targetImg);
        outName = fullfile(pth, 'continuousOverlapMeasures.mat');
        save(outName, 'overlapMeasures');
        
        
       
        % describe output
        aap = aas_desc_outputs(aap, subjind, 'continuousoverlapmeasures', outName);        
end