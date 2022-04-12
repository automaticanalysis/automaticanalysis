function [aap,resp]=aamod_CONN_valid(aap,task)
resp='';

switch task
%     case 'report'
%         localpath = aas_getpath_bydomain(aap,aap.tasklist.currenttask.domain,[subj,sess]);
%         
%         fdiag = dir(fullfile(localpath,'diagnostic_*.jpg'));
%         if isempty(fdiag)
%             streams=aas_getstreams(aap,'output');
%             for streamind=1:length(streams)
%                 % obtain output
%                 outputfnames = aas_getfiles_bystream(aap,aap.tasklist.currenttask.domain,[subj sess],streams{streamind},'output');
%                 
%                 % perform diagnostics
%                 do_diag(outputfnames);
%             end
%             fdiag = dir(fullfile(localpath,'diagnostic_*.jpg'));
%         end
%         
%         for d = 1:numel(fdiag)
%             aap = aas_report_add(aap,subj,'<table><tr><td>');
%             imgpath = fullfile(localpath,fdiag(d).name);
%             aap=aas_report_addimage(aap,subj,imgpath);
%             aap = aas_report_add(aap,subj,'</td></tr></table>');
%         end
    case 'doit'
        % ROI validation - if requested and available
        inp = aas_getstreams(aap,'input');
        useValidROIs = aas_getsetting(aap,'useValidROI');
        if ~isempty(useValidROIs) && useValidROIs && (numel(inp) >= 3) && aas_stream_has_contents(aap,'study',[],inp{3})
            load(aas_getfiles_bystream(aap,'study',[],inp{3}),'ValidROI');
        end
        
        % find common ROI names
        ROIs = cell(1,aas_getN_bydomain(aap,'subject'));
        for subj = 1:aas_getN_bydomain(aap,'subject')
            load(aas_getfiles_bystream(aap,'subject',subj,'settings'),'CONN_x');
            ROIs{subj} = CONN_x.Preproc.variables.names(contains(CONN_x.Preproc.variables.names,'Atlas')); 
            if subj == 1 
                commonROIs = ROIs{subj};
            else
                commonROIs = intersect(commonROIs,ROIs{subj});
            end
        end
                
        % ROI validation
        if exist('ValidROI','var')
            validROIs = cellfun(@(r) sscanf(r,'Atlas.cluster%d'), commonROIs);
            [~,~,indConnectivity] = intersect(ValidROI.ROIval,validROIs);
            commonROIs = commonROIs(indConnectivity);
        end
        
        % determine final indices
        indROIs = cell(1,aas_getN_bydomain(aap,'subject'));
        for subj = 1:aas_getN_bydomain(aap,'subject')
            [~, indROIs{subj}] = intersect(ROIs{subj},commonROIs);
        end
        
        % save valid ROI names
        fnROIs = fullfile(aas_getstudypath(aap),'validROIs.txt');
        writetable(table(commonROIs'),fnROIs,'FileType','text','WriteVariableNames',false);
        aap = aas_desc_outputs(aap,'study',[],'ROInames',fnROIs);
        
        % select CMs' rows and columns of valid ROIs
        for subj = 1:aas_getN_bydomain(aap,'subject')
            for sess = aap.acq_details.selected_sessions
                fnOutSess = cellstr(aas_getfiles_bystream(aap,'session',[subj sess],'connectivity'));
                for f = 1:numel(fnOutSess)
                    load(fnOutSess{f},'sessionCM');
                    sessionCM = sessionCM(indROIs{subj},indROIs{subj});
                    fnOutSess{f} = spm_file(fnOutSess{f},'prefix','valid_');
                    save(fnOutSess{f},'sessionCM');
                end                
                aap = aas_desc_outputs(aap,'session',[subj sess],'connectivity',fnOutSess);
            end
            
            fnOutSubj = cellstr(aas_getfiles_bystream(aap,'subject',subj,'connectivity'));
            for f = 1:numel(fnOutSubj)
                load(fnOutSubj{f},'subjectCM');
                subjectCM = subjectCM(indROIs{subj},indROIs{subj});
                fnOutSubj{f} = spm_file(fnOutSubj{f},'prefix','valid_');
                save(fnOutSubj{f},'subjectCM');
            end
            aap = aas_desc_outputs(aap,'subject',subj,'connectivity',fnOutSubj);
        end
    case 'checkrequirements'
end
