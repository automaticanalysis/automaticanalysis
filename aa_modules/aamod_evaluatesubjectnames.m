% AA initialisation module - evaluate wildcards in subject names
% Performs search using unix ls to convert wildcards into filenames
% Rhodri Cusack MRC CBU Cambridge 2004
% Tibor Auer MRC CBU Cambridge 2012-2013

function [aap,resp]=aamod_evaluatesubjectnames(aap,task,i)

resp='';

switch task
    case 'doit'
        issubj=0;
        if (length(aap.acq_details.subjects)>=i)
            if (~isempty(aap.acq_details.subjects(i).mriname))
                issubj=1;
            end;
            if (strcmp(aap.acq_details.subjects(i).mriname,'missing'))
                aas_log(aap,0,sprintf('MRI from subject number %d is missing, hope you are doing MEG',i));
                return;
            end;
        end;
        if (~issubj)
            aas_log(aap,1,sprintf('No subject name was specified for subject %d\n',i));
        end;
        switch (aap.directory_conventions.remotefilesystem)
            case 'none'
                if (~aas_ismac)
                    cmd=sprintf('cd %s; ls --color=never -d %s',aap.directory_conventions.rawdatadir,aap.acq_details.subjects(i).mriname);
                else
                    cmd=sprintf('cd %s; ls -d %s',aap.directory_conventions.rawdatadir,aap.acq_details.subjects(i).mriname);
                end;
                [s w]=aas_shell(cmd);
                if (s)
                    if isempty(aap.acq_details.subjects(i).megname)
                        aas_log(aap,0,sprintf('Problem finding subject %d raw data directory %s\n',i,aap.acq_details.subjects(i).mriname));
                    else
                        fprintf(' - Warning: Failed to find MRI for subject %d: %s\n',i,aap.acq_details.subjects(i).mriname);
                        aap.acq_details.subjects(i).mriname='missing';
                    end
                else
					% [TA}
                    [a,b] = strtok(w);
                    % on some shells w will have 2 lines (one feedbacking
                    % the new dir after cd, another with the ls result) -
                    % in this case we want the second output from strtok
                    if ~isempty(b) && (numel(b) > 8)
                        a = b;
                    end
                    aap.acq_details.subjects(i).mriname=deblank(strtok(a));
                end;
            case 's3'
                global aaworker
                % Separately match subject and visit parts
                mriname=aap.acq_details.subjects(i).mriname;
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
                    aas_log(aap,true,sprintf('Cannot find raw data for visit matching filter %s. These should now be regular expressions (e.g., CBU090800.*/.*)',aap.acq_details.subjects(i).mriname));
                elseif (length(find_ind)>1)
                    aas_log(aap,true,sprintf('Found more than one raw data set matching visit filter %s. These should now be regular expressions (e.g., CBU090800.*/.*)',aap.acq_details.subjects(i).mriname));
                end;
                aap.acq_details.subjects(i).mriname=s3resp.CommonPrefixes(find_ind).Prefix((length(aap.directory_conventions.rawdatadir)+1):end);
                
                % Check registered in Drupal and put nid in acq_details
                [pth nme ext]=fileparts(aas_getsubjpath(aap,i));
                attr=[];
                attr.datasetalias.value=[nme ext];
                % Check bucket nid
                if (~isfield(aaworker,'bucket_drupalnid'))
                    [aap waserror aaworker.bucket_drupalnid]=drupal_checkexists(aap,'bucket',aaworker.bucket);
                end;
                % Check dataset nid for this subject
                [aap waserror aap.acq_details.subjects(i).drupalnid]=drupal_checkexists(aap,'dataset',aap.acq_details.subjects(i).mriname,attr,aaworker.bucket_drupalnid,aaworker.bucket);

        end;
        
        
    case 'checkrequirements'
        
end;
