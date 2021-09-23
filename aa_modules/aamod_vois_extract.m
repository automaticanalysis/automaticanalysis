% AA module - extract VOIs
% [aap,resp]=aamod_vois_extract(aap,task,subj)
% based on https://jacoblee.net/occamseraser/2018/01/03/extracting-rois-for-ppi-analysis-using-spm-batch/index.html
% Tibor Auer MRC CBSU Jul 2014

function [aap,resp]=aamod_vois_extract(aap,task,subj)

resp='';

switch task
    case 'doit'
        
        %% Init
        job0.spmmat = {aas_getfiles_bystream(aap,'subject',subj,'firstlevel_spm')};
        dat = load(job0.spmmat{1}); SPM = dat.SPM;
        switch aas_getsetting(aap,'adjust')
            case 'none'
                job0.adjust = 0;
            case 'all'
                job0.adjust = nan;
            otherwise
                aas_log(aap,true,'NYI: adjusting for specific contrast')
        end
        imgstream = aas_getstreams(aap,'input'); imgstream = imgstream{end};
        
        for sess = aap.acq_details.selected_sessions   
            aas_makedir(aap,aas_getsesspath(aap,subj,sess));
            
            %% Generate VOIs
            fnames = {};
            for v = aas_getsetting(aap,'VOI')                
                job = job0;
                job.session = sess;
                job.name = fullfile(aas_getsesspath(aap,subj,sess),v.name);
                
                switch v.type
                    case 'sphere'
                        job.roi{1}.sphere.radius = v.size/2;
                    case 'mask'
                        job.roi{1}.mask.image = {aas_getfiles_bystream(aap,'subject',subj,imgstream)};
                        job.roi{1}.mask.threshold = 0.5;
                    case 'roi'
                        job.roi{1}.label.image = {aas_getfiles_bystream(aap,'subject',subj,imgstream)};
                        job.roi{1}.label.list = v.centredefinition.roival;
                    case 'blob'
                        aas_log(aap,true,'NYI: blob')
                end
                switch v.centre
                    case 'xyz'
                        aas_log(aap,true,'NYI: xyz')
                    case 'roicentre'
                        [Y,XYZmm] = spm_read_vols(spm_vol(aas_getfiles_bystream(aap,'subject',subj,imgstream)));
                        job.roi{1}.(char(fieldnames(job.roi{1}))).centre = mean(XYZmm(:,Y(:) == v.centredefinition.roival),2)';
                    case 'roimaximum'
                        job.roi{1}.(char(fieldnames(job.roi{1}))).centre = [0 0 0]; % not used
                        job.roi{1}.(char(fieldnames(job.roi{1}))).move.local.spm = 2;
                        job.roi{1}.(char(fieldnames(job.roi{1}))).move.local.mask = 'i3';
                        job.roi{2}.spm.spmmat = {''};
                        job.roi{2}.spm.contrast = find(strcmp({SPM.xCon.name},v.centredefinition.contrast),1,'first');
                        job.roi{2}.spm.conjunction = 1;
                        job.roi{2}.spm.mask = [];
                        job.roi{2}.spm.threshdesc = 'none';
                        job.roi{2}.spm.thresh = 1; % unthresholded
                        job.roi{2}.spm.extent = 0;
                        job.roi{3}.label.image = {aas_getfiles_bystream(aap,'subject',subj,imgstream)};
                        job.roi{3}.label.list = v.centredefinition.roival;
                        job.expression = 'i1';
                end
                
                % add brainmask
                job.roi{end+1}.mask.image = {aas_getfiles_bystream(aap,'subject',subj,'firstlevel_brainmask')};
                job.roi{end}.mask.threshold = 0.5;
                job.expression = sprintf('%s & i%d',job.expression,numel(job.roi));
                
                out = spm_run_voi(job);
                fnames = vertcat(fnames,out.voiimg,out.voimat);
            end
            aap = aas_desc_outputs(aap,'session',[subj sess],'vois',fnames);
        end
        
        
    case 'checkrequirements'
end
