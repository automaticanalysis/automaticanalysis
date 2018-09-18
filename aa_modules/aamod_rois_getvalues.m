% AA module - UNNORMALISE_ROI Unnormalise a particular ROI for a given subject
%
% Alejandro Vicente-Grabovetsky Dec-2008

function [aap,resp] = aamod_rois_getvalues( aap, task, subj, sess)

resp='';

switch task
    case 'doit'
        
        ROIlist = aas_getfiles_bystream(aap, subj, 'rois');
        
        streamsIn = aap.tasklist.currenttask.inputstreams.stream;
        streamsIn = streamsIn(~strcmp(streamsIn, 'rois'));
        
        if length(streamsIn) > 1
           aas_log(aap, true, 'aamod_rois_getvalues can only handle 1 stream of data') 
        end
        
        % Is session specified in task header?
        if (isfield(aap.tasklist.currenttask.settings,'session'))
            sess = aap.tasklist.currenttask.settings.session;
            dataList = aas_getfiles_bystream(aap, subj, sess, streamsIn{1});
        else
            dataList = aas_getfiles_bystream(aap, subj, streamsIn{1});
        end
        
        for f = size(dataList):-1:1
            [pth, nme, ext] = fileparts(deblank(dataList(f,:)));
            if strcmp(ext, '.hdr')
                dataList(f,:) = [];
            end
        end
        
        roiValues = nan(size(ROIlist, 1), size(dataList, 1));
        
        % Loop through all ROIs
        for r = 1:size(ROIlist, 1)
            rV = spm_vol(deblank(ROIlist(r,:)));
            rY = spm_read_vols(rV);
            
            for f = 1:size(dataList)
                dV = spm_vol(deblank(dataList(f,:)));
                if any(dV.dim ~= rV.dim)
                    aas_log(aap, true, 'Your ROI matrix size does not correspond to your data matrix size')
                end
                dY = spm_read_vols(dV);
                
                 vals = dY(logical(rY));
                 vals = vals(~isnan(vals));
                 switch aap.tasklist.currenttask.settings.metric
                     case 'mean'
                        roiValues(r, f) = mean(vals);
                     case 'median'
                        roiValues(r, f) = median(vals);
                     case 'trimmean'
                         roiValues(r, f) = trimmean(vals, 95);
                     case 'robustfit'
                         roiValues(r, f) = robustfit(zeros(length(vals), 0), vals);                         
                 end
            end
        end
        
        streamsOut = aap.tasklist.currenttask.outputstreams.stream;
        if (isfield(aap.tasklist.currenttask.settings,'session'))
            outstream = fullfile(aas_getsesspath(aap, subj, sess), 'roiValues.mat');
            save(outstream, 'roiValues');
            aap=aas_desc_outputs(aap, subj, sess, streamsOut{1}, outstream);
        else
            outstream = fullfile(aas_getsubjpath(aap, subj), 'roiValues.mat');
            save(outstream, 'roiValues');
            aap = aas_desc_outputs(aap, subj, streamsOut{1}, outstream);
        end
end