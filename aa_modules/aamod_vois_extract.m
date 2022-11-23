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
        
        contrast_specifier = aas_getsetting(aap,'adjust');
        
        switch contrast_specifier
            
            case 'none'
                
                job0.adjust = 0;
                
            case 'all'
                
                job0.adjust = nan;
                
            otherwise
                
                % adjust by a named contrast (e.g., an F "effects of interest" contrast)
                                     
                job0.adjust = find(strcmp({SPM.xCon.name},contrast_specifier),1,'first'); 
                
                if isempty(job0.adjust)
                    aas_log(aap,true,sprintf('Invalid <adjust> parameter %s',contrast_specifier));                      
                end
                 
        end

        % the rois instream (if present) can either be used as a
        % mask (type=mask) or an atlas (type=roi). 

        rois_instream = [];
        if (aas_stream_has_contents(aap, subj, 'rois'))
            rois_instream = aas_getstreams(aap,'input'); rois_instream = rois_instream{end}; 
        end
        
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
                        
                        if isempty(v.size)
                            aas_log(aap,true,'You must specify a diameter in <size> to use the type=sphere option');                      
                        end
                     
                        job.roi{1}.sphere.radius = v.size/2;

                    case 'mask'

                        if isempty(rois_instream)       
                            aas_log(aap,true,'You must specify a mask in the rois instream to use the type=mask option');
                        end
                        
                        if isempty(v.maskthresh)       
                            aas_log(aap,true,'You must specify a mask threshold in <maskthresh> to use the type=mask option');
                        end
                      
                        job.roi{1}.mask.image = {aas_getfiles_bystream(aap,'subject', subj, rois_instream)};
                        job.roi{1}.mask.threshold = v.maskthresh;                    

                    case 'roi'
                        
                        if isempty(rois_instream)       
                            aas_log(aap,true,'You must specify an atlas in the rois instream to use the type=roi option');
                        end
                        
                        if isempty(v.centredefinition.roival)
                            aas_log(aap,true,'You must specify an atlas label in <centredefintion><roival> to use the type=roi option');                      
                        end
                        
                        job.roi{1}.label.image = {aas_getfiles_bystream(aap,'subject', subj, rois_instream)};
                        job.roi{1}.label.list = v.centredefinition.roival;
                                                
                end
                
                if ~isempty(v.centre) % mask-type ROI definitions don't include a centre

                    if ~strcmp(v.type,'sphere')
                        aas_log(aap,true,'You must use type=sphere when specifying a <centre> option');                      
                    end

                    switch v.centre

                        case 'xyz'

                            % we interpret center=xyz specifies a fixed sphere -- (v.type = sphere assumed)

                            if isempty(v.centredefinition.xyz)
                                aas_log(aap,true,'You must specify a sphere center in <centredefintion><xyz> to use the centre=xyz option');                      
                            end

                            job.roi{1}.sphere.centre = v.centredefinition.xyz;
                            job.roi{1}.sphere.move.fixed = 1;

                        case 'roicentre'

                            if isempty(rois_instream)       
                                aas_log(aap,true,'You must specify an atlas in the rois instream to use the centre=roicentre option');
                            end

                            if isempty(v.centredefinition.roival)
                                aas_log(aap,true,'You must specify an atlas label in <centredefintion><roival> to use the centre=roicentre option');                      
                            end

                            [Y,XYZmm] = spm_read_vols(spm_vol(aas_getfiles_bystream(aap,'subject',subj,rois_instream)));
                            job.roi{1}.(char(fieldnames(job.roi{1}))).centre = mean(XYZmm(:,Y(:) == v.centredefinition.roival),2)';

                            % this option assumes sphere is fixed

                            job.roi{1}.sphere.move.fixed = 1; 


                        case 'roimaximum'

                            if isempty(rois_instream)       
                                aas_log(aap,true,'You must specify an atlas in the rois instream to use the centre=roimaximum option');
                            end

                            if isempty(v.centredefinition.roival)
                                aas_log(aap,true,'You must specify an atlas label in <centredefintion><roival> to use the centre=roimaximum option');                      
                            end

                             if isempty(v.centredefinition.contrast)
                                aas_log(aap,true,'You must specify a contrast in <centredefintion><contrast> to use the centre=roimaximum option');                      
                             end

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
                            job.roi{3}.label.image = {aas_getfiles_bystream(aap,'subject',subj,rois_instream)};
                            job.roi{3}.label.list = v.centredefinition.roival;

                    end
                    
                end
                
                % we always mask everything against the firstlevel brainmask
                              
                job.roi{end+1}.mask.image = {aas_getfiles_bystream(aap,'subject',subj,'firstlevel_brainmask')};
                job.roi{end}.mask.threshold = 0.5; % hardcode an appropriate threshold for a 0/1 mask
              
                % note you MUST pass an expression or spm_run_voi will crash
                job.expression = sprintf('i1 & i%d',numel(job.roi));
                
                out = spm_run_voi(job);
                
                fnames = vertcat(fnames,out.voiimg,out.voimat);
                
            end
            aap = aas_desc_outputs(aap,'session',[subj sess],'vois',fnames);
        end
        
        
    case 'checkrequirements'
end
