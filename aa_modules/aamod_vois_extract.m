% AA module -
% function [aap,resp]=aamod_vois_extract(aap,task,i,j)
% Rhodri Cusack 2006-2012

function [aap,resp]=aamod_vois_extract(aap,task,i,j)

resp='';

switch task
    case 'doit'
        vois=aap.tasklist.currenttask.settings.vois;
        for r=1:length(vois)
            YROI=[];
            % which files?
            files = aas_getimages_bystream(aap,i,j,'epi');
            
            % Define ROI if not already done
            if (isempty(YROI))
                Vepi=spm_vol(files(1,:));
                switch vois(r).type
                    case 'sphere'
                        [Yepi XYZ]=spm_read_vols(Vepi);
                        XYZ=XYZ';
                        YROI=false(size(Yepi));
                        roiradsq=vois(r).radius ^ 2;
                        for reg=1:size(vois(r).centres,1)
                            d=XYZ-repmat(vois(r).centres(reg,:),[size(XYZ,1) 1]);
                            d=sum(d.^2,2);
                            YROI=YROI | reshape(d<=roiradsq,size(Yepi));
                        end;
                    case 'image'
                        VROI=spm_vol(vois(r).filename);
                        if (any(VROI.dim~=Vepi.dim))
                            aas_log(aap,1,fprintf('VOI images must have same dimensions (i.e., matrix size) as EPIs. Not true for VOI %s\n',vois(r).filename));
                        end;
                        if (any(VROI.mat~=Vepi.mat))
                            aas_log(aap,1,fprintf('VOI images must have space (voxel sizes, centre in .mat file) as EPIs. Not true for VOI %s\n',vois(r).filename));
                        end;
                        YROI=spm_read_vols(VROI);
                    otherwise
                        aas_log(aap,1,fprintf('VOI type %s unknown (VOI %s)\n',vois(r).type,vois(r).name));
                end;
            end;
            
            [x y z]=ind2sub(size(YROI),find(YROI>0));
            XYZROI=[x,y,z];
            % now extract data
            y=spm_get_data(files,XYZROI');
        
            % compute mean response
            av=mean(y');
            
            % compute regional response in terms of first eigenvariate
            nev=aap.tasklist.currenttask.settings.numberofeigenvalues;
            [m n]   = size(y);
            if m > n
                [v s v] = svd(y'*y);
                s       = diag(s);
                v       = v(:,1:nev);
                u       = y*v./repmat(sqrt(s(1:nev)'),[size(y,2) 1]);
            else
                [u s u] = svd(y*y');
                s       = diag(s);
                u       = u(:,1:nev);
                v       = y'*u./repmat(sqrt(s(1:nev)'),[size(y,2) 1]);
            end
            d       = sign(sum(v));
            u       = u.*repmat(d,[size(u,1) 1]);
            v       = v.*repmat(d,[size(v,1) 1]);
            Y       = u.*repmat(sqrt(s(1:nev)'/n),[size(u,1) 1]);
            
            vois(r).data.raw=y;
            if (nev==1)
                vois(r).data.firsteigenvalue=Y;
            else
                vois(r).data.topeigenvalues=Y;
            end;
            vois(r).data.mean=y;
            
            % Need to understand how this differs from the eigenvariate
            % above [RC]
            [coeff score latent]=princomp(vois(r).data.raw);
            vois(r).pca.coeff=coeff;
            vois(r).pca.score=score;
            vois(r).pca.latent=latent;
            vois(r).xyz=XYZROI;
            
        end;
        
        % and save
        outdir=aas_getsesspath(aap,i,j);
        aas_makedir(aap,outdir);
        
        outfn=fullfile(outdir,sprintf('vois_%s.mat',strtrim(aap.tasklist.currenttask.settings.name)));
        vois=vois;
        save(outfn,'vois');
        aap=aas_desc_outputs(aap,i,j','voi',outfn);
        
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end;