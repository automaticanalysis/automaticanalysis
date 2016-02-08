% AA initialisation module - evaluate wildcards in subject names
% Performs search using unix ls to convert wildcards into filenames
% Rhodri Cusack MRC CBU Cambridge 2004


function [aap,resp]=aamod_evaluatesubjectnames(aap,task,subj)

resp='';

switch task
    case 'doit'
        issubj=false;
        isMRI = false; isMEG = false;
        if ~iscell(aap.acq_details.subjects(subj).mriname), aap.acq_details.subjects(subj).mriname = {aap.acq_details.subjects(subj).mriname}; end
        
        if numel(aap.acq_details.subjects)>=subj
            if ~isempty(cell2mat(aap.acq_details.subjects(subj).mriname))
                issubj = true;
                isMRI = true;
            end;
            if all(strcmp(aap.acq_details.subjects(subj).mriname,'missing'))
                aas_log(aap,0,sprintf('MRI from subject number %d is missing, hope you are doing MEG',subj));
                issubj = false;
                isMRI = false;
            end;
            if ~isempty(aap.acq_details.subjects(subj).megname)
                issubj = true;
                isMEG = true;
            end
        end;
        if ~issubj
            aas_log(aap,true,sprintf('No subject name was specified for subject %d\n',subj));
        end;
        switch (aap.directory_conventions.remotefilesystem)
            case 'none'
                if isMRI
                    for m = 1:numel(aap.acq_details.subjects(subj).mriname)
                        s = mri_findvol(aap,aap.acq_details.subjects(subj).mriname{m});
                        if isempty(s)
                            aas_log(aap,false,sprintf('WARNING: Problem finding raw MRI data directory for subject %s',aap.acq_details.subjects(subj).subjname));
                            if ~isempty(aap.acq_details.subjects(subj).megname)
                                aas_log(aap,false,'INFO: You can still run MEG analysis');
                                aap.acq_details.subjects(subj).mriname{m}='missing';
                            end
                        else
                            aap.acq_details.subjects(subj).mriname{m}=s;
                        end;
                    end
                    iname = find(~strcmp(aap.acq_details.subjects(subj).mriname,'missing'),1);
                    if ~isempty(iname)
                        name = aas_mriname2subjname(aap,aap.acq_details.subjects(subj).mriname{iname});
                    end
                end
                if isMEG
                    s = meg_findvol(aap,aap.acq_details.subjects(subj).megname);
                    if isempty(s)
                        if isempty(aap.acq_details.subjects(subj).mriname) || all(strcmp(aap.acq_details.subjects(subj).mriname,'missing'))
                            aas_log(aap,true,sprintf('ERROR: Problem finding raw MEG data directory for subject %s',aap.acq_details.subjects(subj).subjname));
                        else
                            aas_log(aap,false,'WARNING: Problem finding raw MEG data directory for subject %s\n',aap.acq_details.subjects(subj).subjname);
                            aas_log(aap,false,'INFO: You can still run MRI analysis');
                            aap.acq_details.subjects(subj).megname='missing';
                        end
                    else
                        aap.acq_details.subjects(subj).megname=s;                        
                    end
                    iname = find(~strcmp(aap.acq_details.subjects(subj).megname,'missing'),1);
                    if ~isempty(iname)
                        name = aas_megname2subjname(aap,aap.acq_details.subjects(subj).megname);
                    end
                end
                if ~isfield(aap.acq_details.subjects(subj),'subjname') || isempty(aap.acq_details.subjects(subj).subjname)
                    switch aap.directory_conventions.subject_directory_format
                        case 0 % from predefined list
                            aap.acq_details.subjects(subj).subjname = aap.directory_conventions.subject_directory_names{subj};
                        case 1 % from data
                            aap.acq_details.subjects(subj).subjname = name;
                        case 2 % S#
                            aap.acq_details.subjects(subj).subjname = sprintf('S%02d',subj);
                        case 3 % full path/manual
                            aap.acq_details.subjects(subj).subjname = aap.directory_conventions.subject_directory_names{subj};
                        otherwise
                            aas_log(aap,true,sprintf('Unknown subject directory format (aap.directory_conventions.subject_directory_format=%d',aap.directory_conventions.subject_directory_format));
                    end
                end
            case 's3' 
                % [TA] needs changing to enable 
                %   - multiple rawdatadir
                %   - multiple mriname
                %   - use of subjname
                global aaworker
                % Separately match subject and visit parts
                mriname=aap.acq_details.subjects(subj).mriname;
                while (mriname(end)=='/')
                    mriname=mriname(1:(end-1));
                end;
                [pth nme ext]=fileparts(mriname);
                subjectfilter=pth;
                visitfilter=[nme ext];
                
                % First subject, get list of all subjects
                [aap s3resp]=s3_list_objects(aap,aaworker.bucketfordicom,[aap.directory_conventions.rawdatadir '/'],[],'/');
                ind=cellfun(@(x) ~isempty(regexp(x,fullfile(aap.directory_conventions.rawdatadir,subjectfilter))),{s3resp.CommonPrefixes.Prefix});
                find_ind=find(ind);
                if (isempty(find_ind))
                    aas_log(aap,true,sprintf('Cannot find raw data for subject matching filter %s. These should now be regular expressions (e.g., CBU090800.*/.*)',subjectfilter));
                elseif (length(find_ind)>1)
                    aas_log(aap,true,sprintf('Found more than one raw data set matching subject filter %s. These should now be regular expressions (e.g., CBU090800.*/.*)',subjectfilter));
                end;
                matching_subject=s3resp.CommonPrefixes(find_ind).Prefix; % this will be the full path with trailing /... note this below
                
                % Now, visit filter, get list of all visits
                [aap s3resp]=s3_list_objects(aap,aaworker.bucketfordicom,matching_subject,[],'/');
                ind=cellfun(@(x) ~isempty(regexp(x,visitfilter)),{s3resp.CommonPrefixes.Prefix});
                find_ind=find(ind);
                if (isempty(find_ind))
                    aas_log(aap,true,sprintf('Cannot find raw data for visit matching filter %s. These should now be regular expressions (e.g., CBU090800.*/.*)',aap.acq_details.subjects(subj).mriname));
                elseif (length(find_ind)>1)
                    aas_log(aap,true,sprintf('Found more than one raw data set matching visit filter %s. These should now be regular expressions (e.g., CBU090800.*/.*)',aap.acq_details.subjects(subj).mriname));
                end;
                aap.acq_details.subjects(subj).mriname=s3resp.CommonPrefixes(find_ind).Prefix((length(aap.directory_conventions.rawdatadir)+1):end);
                
                % Check registered in Drupal and put nid in acq_details
                [pth nme ext]=fileparts(aas_getsubjpath(aap,subj));
                attr=[];
                attr.datasetalias.value=[nme ext];
                % Check bucket nid
                if (~isfield(aaworker,'bucket_drupalnid'))
                    [aap waserror aaworker.bucket_drupalnid]=drupal_checkexists(aap,'bucket',aaworker.bucket);
                end;
                % Check dataset nid for this subject
                [aap waserror aap.acq_details.subjects(subj).drupalnid]=drupal_checkexists(aap,'dataset',aap.acq_details.subjects(subj).mriname,attr,aaworker.bucket_drupalnid,aaworker.bucket);

        end;
                
    case 'checkrequirements'
        
end;
