function [aap,resp] = aamod_roi_valid(aap,task)

% AAMOD_ROI_VALID Look for valid ROIs (with non-NaN values across subjects)
%
% INPUT options [defaults if not set in xml file]:
%  inputstreams
%   .stream{1}    = ROI data stream for input [roidata_epi]
%  AbsVoxThr      = Threshold min N voxels [10]
%  SubjRemoveStat = Central tendency statistic for subject threshold [mode]
%  SubjRemoveThr  = Variability statistic multiplier for subject threshold []
%  ROIRemoveThr   = Threshold (subject count) for ROI removal [0]
%
% OUTPUT:
%  outputstreams
%   .stream{1}   = ROI valid list [roivalid]
%
% by Jason Taylor (09/Aug/2013)

resp='';

switch task
    
    case 'domain'
        resp='study';
        
    case 'description'
        resp='Find valid ROIs across subjects';
        
    case 'doit'
        
        try AbsVoxThr      = aap.tasklist.currenttask.settings.AbsVoxThr;      catch, AbsVoxThr      = 10; end
        try SubjRemoveStat = aap.tasklist.currenttask.settings.SubjRemoveStat; catch, SubjRemoveStat = 0.8; end
        try SubjRemoveThr  = aap.tasklist.currenttask.settings.SubjRemoveThr;  catch, SubjRemoveThr  = 3;  end
        try ROIRemoveThr   = aap.tasklist.currenttask.settings.ROIRemoveThr;   catch, ROIremoveThr   = 0;  end
        
        instreams = aas_getstreams(aap,'input');
        for i = 1:numel(instreams)
            instream = instreams{i};
            
            %% Do it:
            ValidROI = struct();
            Nv       = [];
            Mm       = [];
            ROIval   = [];
            
            sourcedomain = aap.internal.inputstreamsources{aap.tasklist.currenttask.modulenumber}.stream(i).sourcedomain;
            if strcmp(sourcedomain,'subject')
                procind = 1;
                indind = 1;
            elseif ~isempty(strfind(sourcedomain,'session'))
                procind = aap.acq_details.selected_sessions;
                indind = 1:2;
            end
            
            for p = procind

                for subjind = 1:length(aap.acq_details.subjects),
                    
                    % Load ROI file for subject/session:
                    indices = [subjind procind];
                    ROIfname = aas_getfiles_bystream(aap,sourcedomain,indices(indind),instream);
                    loaded = load(ROIfname); ROI = loaded.ROI;
                    % Get number of valid voxels in each ROI:
                    Nv(subjind,:) = [ROI.Nvox];
                    mROI = [ROI.mean];
                    Mm(subjind,:) = mROI;
                    
                end
                
                [pth stem fext] = fileparts(ROIfname);
                invalidroi = isnan(Mm) | Nv<AbsVoxThr;
                ROIval = [ROI.ROIval];
                Nr = length(ROIval);
                
                
                % Remove bad subjects before removing ROIs:
                
                switch SubjRemoveStat,
                    case 'mode'
                        % Keep only if N(invalid) is mode or less:
                        scrit = mode(sum(invalidroi,2));
                    case 'median'
                        % Keep only if N(invalid) is median+Thr*IQR or less:
                        scrit = median(sum(invalidroi,2));
                        scrit = scrit + SubjRemoveThr*iqr(sum(invalidroi,2));
                    case 'mean'
                        % Keep only if N(invalid) is mean+Thr*STD or less:
                        scrit = mean(sum(invalidroi,2));
                        scrit = scrit + SubjRemoveThr*std(sum(invalidroi,2));
                end
                
                % Find rubbish subjects:
                s2ignore = sum(invalidroi,2)>scrit;
                
                % Find ROIs to ignore (ignoring rubbish subjects):
                r2ignore = sum(invalidroi(~s2ignore,:))>ROIRemoveThr;
                
                % Session Summary:
                ValidROI(p).sessname          = aas_getsessname(aap,p);
                ValidROI(p).ROIval            = ROIval(~r2ignore);
                ValidROI(p).ROIind            = setdiff(1:length(ROI),find(r2ignore));
                ValidROI(p).ROIind2ignore     = find(r2ignore);
                ValidROI(p).Subjind           = setdiff(1:length(aap.acq_details.subjects),find(s2ignore));
                ValidROI(p).Subjind2ignore    = find(s2ignore);
                ValidROI(p).AbsVoxThr         = AbsVoxThr;
                ValidROI(p).SubjRemoveStat    = SubjRemoveStat;
                ValidROI(p).SubjRemoveThr     = SubjRemoveThr;
                ValidROI(p).SubjRemoveCritVal = scrit;
                ValidROI(p).ROIRemoveThr      = ROIRemoveThr;
                
                % plot imagesc:
                %f=spm_figure('FindWin'); spm_figure('Clear',f);
                %imagesc(invalidroi); colormap gray;
                [pth,nm] = fileparts(ROIfname);
                ind=strfind(pth,'/');
                ind=ind(end)-1; % remove subject, session
                pth = pth(1:ind);
                %outfile = fullfile(pth,sprintf('imagesc_ValidROI_%s.png',aap.acq_details.sessions(sessind).name));
                %print(f,'-dpng',outfile);
                %aap = aas_desc_outputs(aap,sprintf('imagesc_%s_%s',aap.acq_details.sessions(sessind).name,outstream),outfile);
                
            end
            
            %% Describe output
            outstream = strrep(instream,'roidata','valid_roi');
            outfile = fullfile(aas_getstudypath(aap),[strrep(outstream,'valid_roi','ValidROI') '.mat']);
            save(outfile,'ValidROI');
            aap = aas_desc_outputs(aap,'study',[],outstream,outfile);
        end
        
    case 'checkrequirements'
        % get input
        [stagename, index] = strtok_ptrn(aap.tasklist.currenttask.name,'_0');
        stageindex = sscanf(index,'_%05d');
        in = aap.tasksettings.(stagename)(stageindex).inputstreams.stream; if ~iscell(in), in = {in}; end
        
        % switch for source stage
        aaps = aas_setcurrenttask(aap,aap.internal.inputstreamsources{aap.tasklist.currenttask.modulenumber}.stream.sourcenumber);
        src = aas_getstreams(aaps,'output');
        for s = numel(in)
            if ~strcmp(in{s},src{s})
                aap = aas_renamestream(aap,aap.tasklist.currenttask.name,in{s},src{s},'input');
            end
            aas_log(aap,false,['INFO: ' aap.tasklist.currenttask.name ' input stream: ''' src{s} '''']);
        end
end
