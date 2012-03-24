% AA module - first level statistics
% function [aap,resp]=aamod_vois_extract(aap,task,i)
% You will need to make a local copy of this module and then modify it to
%  reflect the particular experiment you conducted
% First-level model based on original by Adam Hampshire MRC CBU Cambridge Feb 2006
% Modified for aa by Rhodri Cusack MRC CBU Mar-May 2006
% @@@ THIS IS NOT YET TRANSFORMED TO AA4 @@@

function [aap,resp]=aamod_vois_extract(aap,task,i)

resp='';

switch task
    case 'domain'
        resp='subject';   % this module needs to be run once per subject

    case 'description'
        resp='SPM5 VOI extract';

    case 'summary'
        subjpath=aas_getsubjpath(i);
        resp=sprintf('SPM5 VOI extract for %s\n',subjpath);

    case 'report'

    case 'doit'

        for r=1:length(aap.spmanalysis.vois)
            YROI=[];
            for sess = aap.acq_details.selected_sessions
                % which files?
                files = aas_getimages(aap,i,sess,aap.tasklist.currenttask.epiprefix);                
                
                % Define ROI if not already done
                if (isempty(YROI))
                    Vepi=spm_vol(files(1,:));
                    switch aap.spmanalysis.vois{r}.type
                        case 'sphere'
                            [Yepi XYZ]=spm_read_vols(Vepi);
                            XYZ=XYZ';
                            YROI=false(size(Yepi));
                            roiradsq=aap.spmanalysis.vois{r}.radius ^ 2;
                            for reg=1:size(aap.spmanalysis.vois{r}.centres,1)
                                d=XYZ-repmat(aap.spmanalysis.vois{r}.centres(reg,:),[size(XYZ,1) 1]);
                                d=sum(d.^2,2);
                                YROI=YROI | reshape(d<=roiradsq,size(Yepi));
                            end;
                        case 'image'
                            VROI=spm_vol(aap.spmanalysis.vois{r}.filename);
                            if (any(VROI.dim~=Vepi.dim))
                                aas_log(aap,1,fprintf('VOI images must have same dimensions (i.e., matrix size) as EPIs. Not true for VOI %s\n',aap.spmanalysis.vois{r}.filename));
                            end;
                            if (any(VROI.mat~=Vepi.mat))
                                aas_log(aap,1,fprintf('VOI images must have space (voxel sizes, centre in .mat file) as EPIs. Not true for VOI %s\n',aap.spmanalysis.vois{r}.filename));
                            end;
                            YROI=spm_read_vols(VROI);
                        otherwise
                            aas_log(aap,1,fprintf('VOI type %s unknown (VOI %s)\n',aap.spmanalysis.vois{r}.type,aap.spmanalysis.vois{r}.name));
                    end;
                end;
                
                [x y z]=ind2sub(size(YROI),find(YROI>0));
                XYZROI=[x,y,z];
                % now extract data
                y=spm_get_data(files,XYZROI');
                    
                % compute mean response
                av=mean(y');
                
                % compute regional response in terms of first eigenvariate
                [m n]   = size(y);
                if m > n
                    [v s v] = svd(spm_atranspa(y));
                    s       = diag(s);
                    v       = v(:,1);
                    u       = y*v/sqrt(s(1));
                else
                    [u s u] = svd(spm_atranspa(y'));
                    s       = diag(s);
                    u       = u(:,1);
                    v       = y'*u/sqrt(s(1));
                end
                d       = sign(sum(v));
                u       = u*d;
                v       = v*d;
                Y       = u*sqrt(s(1)/n);
                
                aap.spmanalysis.vois{r}.data{sess}.raw=y;
                aap.spmanalysis.vois{r}.data{sess}.firsteigenvalue=Y;
                aap.spmanalysis.vois{r}.data{sess}.mean=y;

            end;

            % and save
            outdir=fullfile(aas_getsubjpath(aap,i),aap.directory_conventions.voidirname);
            aas_makedir(aap,outdir);
            outfn=fullfile(outdir,aap.directory_conventions.voifilename);
            vois=aap.spmanalysis.vois;  
            save(outfn,'vois');
                
        end;

    case 'checkrequirements'

    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end;