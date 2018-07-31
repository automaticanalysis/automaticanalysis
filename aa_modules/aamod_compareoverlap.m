function [aap,resp] = aamod_compareoverlap(aap, task, subjind)
% AAMOD_COMPAREOVERLAP Compare overlap of two binary images
%
% [aap,resp] = AAMOD_COMPAREOVERLAP(aap, task, subjind)
%
% input streams:   2 images (first "target", second "test")
% output stream:   overlapmeasures
%
% Target image assumed to be binary. Input image assumed continuous but can
% be binary:
%
%  aap.tasklist.currenttask.settings.binary = 1;
%
% Both are assumed to be in the same space and voxel-by-voxel matched.
%
% If thresh specified, binarize test image at one or more thresholds. E.g.,
%
%  aap.tasklist.currenttask.settings.thresh = [.7 .8 .9];
%
% will produce three sets of comparisons for different thresholds.
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



resp='';

switch task
    case 'domain'
        resp='subject';
    case 'description'
        resp='Compare overlap of test image with a target image';
    case 'doit'
                        
        targetStream = aap.tasklist.currenttask.inputstreams(1).stream{1};
        testStream = aap.tasklist.currenttask.inputstreams(1).stream{2};        
        
        targetImg = aas_getfiles_bystream(aap, subjind, targetStream);
        testImg = aas_getfiles_bystream(aap, subjind, testStream);
        
        testIsBinary = aap.tasklist.currenttask.settings.binary;
        
        % Threshold(s) for binarizing. To disable, set to Inf or [].
        thresh = aap.tasklist.currenttask.settings.thresh;
        if isempty(thresh) || isinf(thresh(1))
            thresholdImage = false;
        else
            thresholdImage = true;
        end
        
        
        % If test image is continuous, compute continuous overlaps. If
        % binary, don't bother.
        if testIsBinary
            computeContinuous = false;
        else
            computeContinuous = true;
        end
        
        % If input image is binary, or we are binarizing using
        % thresholding, compute binary measures.
        if testIsBinary || thresholdImage
            computeBinary = true;
        else
            computeBinary = false;
        end
        
        
        % read in images and reshape to be 1 by n vector
        Vtarget = spm_vol(targetImg);
        [Ytarget, XYZtarget] = spm_read_vols(Vtarget);
        Ytarget = reshape(Ytarget, 1, numel(Ytarget));
        targetTrueInd = find(Ytarget==1);
        targetFalseInd = find(Ytarget==0); 
        
        Vtest = spm_vol(testImg);
        [Ytest, XYZtest] = spm_read_vols(Vtest);
        Ytest = reshape(Ytest, 1, numel(Ytest));
        
        
        % structure for saving overlap measures
        overlapMeasures = [];
        overlapMeasures.thresh = thresh;
        
        
        if computeContinuous 
            
            % Target overlap
            overlapMeasures.targetOverlapContinuous = sum(Ytest.*Ytarget)./sum(Ytarget);
            
            % Mean overlap / Dice coefficient
            overlapMeasures.meanOverlapContinuous = 2 * (sum(Ytest.*Ytarget)./sum(Ytest + Ytarget));
            
            % Union overlap / Jaccard coefficient
            overlapMeasures.unionOverlapContinuous = (sum(Ytest.*Ytarget)./sum(Ytest + Ytarget));
            
            imageDiff = Ytest-Ytarget;
            
            % MSE
            overlapMeasures.mseContinuous = mean(imageDiff.^2);
            
            % False negative (where imageDiff is negative)
            whereNegative = find(imageDiff<0);
            overlapMeasures.falseNegativeContinuous = sum(Ytest(whereNegative)); % 0 in test, 1 in target
            
            % False positive
            wherePositive = find(imageDiff>0);
            overlapMeasures.falsePostiveContinuous = sum(Ytest(wherePositive)); % 1 in test, 0 in target
            
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
            
            overlapMeasures.dprimeContinuous = norminv(hits, 0, 1) - norminv(falseAlarms, 0, 1);
        end
        
        
        if computeBinary
            if thresholdImage
                for threshInd = 1:length(thresh)
                    YtestBinary = Ytest > thresh(threshInd);
                                        
                    % Target overlap
                    overlapMeasures.targetOverlapBinary(threshInd) = sum(YtestBinary.*Ytarget)./sum(Ytarget);
                    
                    % Mean overlap / Dice coefficient
                    overlapMeasures.meanOverlapBinary(threshInd) = 2 * (sum(YtestBinary.*Ytarget)./sum(YtestBinary + Ytarget));
                    
                    % Union overlap / Jaccard coefficient
                    overlapMeasures.unionOverlapBinary(threshInd) = (sum(YtestBinary.*Ytarget)./sum(YtestBinary + Ytarget));
                    
                    % False negative
                    overlapMeasures.falseNegativeBinary(threshInd) = sum(YtestBinary < Ytarget); % 0 in test, 1 in target
                    
                    % False positive
                    overlapMeasures.falsePostiveBinary(threshInd) = sum(YtestBinary > Ytarget); % 1 in test, 0 in target
                end
            else
                % (No need to threshold, already binary)
                
                % Target overlap
                overlapMeasures.targetOverlapBinary = sum(Ytest.*Ytarget)./sum(Ytarget);
                
                % Mean overlap / Dice coefficient
                overlapMeasures.meanOverlapBinary = 2 * (sum(Ytest.*Ytarget)./sum(Ytest + Ytarget));
                
                % Union overlap / Jaccard coefficient
                overlapMeasures.unionOverlapBinary = (sum(Ytest.*Ytarget)./sum(Ytest + Ytarget));
                
                % False negative
                overlapMeasures.falseNegativeBinary = sum(Ytest < Ytarget); % 0 in test, 1 in target
                
                % False positive
                overlapMeasures.falsePostiveBinary = sum(Ytest > Ytarget); % 1 in test, 0 in target
                                
            end % checking for thresholding
        end
        
        % get filepath and save overlapMeasures.mat                
        [pth, nm, ext] = fileparts(testImg);
        outName = fullfile(pth, 'overlap.mat');
        save(outName, 'overlapMeasures');
        
               
        % describe output
        aap = aas_desc_outputs(aap, subjind, 'overlapmeasures', outName);        
end