% AA module
% Performs DICOM to NIFTI conversion of all data
% You will need modified version of spm_dicom_convert.m available from aa
% website
% Also moves dummies to 'dummy_scans' directory within each session
% Rhodri Cusack MRC CBU Cambridge Nov 2005

function [aap,resp]=aamod_convert_epis(aap,task,i,j)

resp='';

switch task
    case 'description'
        resp=sprintf('Performing dicom to nifti conversion of EPIs');
        
    case 'summary'
        resp=sprintf('Made directory %s and did dicom to nifti conversion of EPIs\n');
        
    case 'report'
        
    case 'doit'
        subjpath=aas_getsubjpath(aap,i);
        % If doesn't exist then make session directory, as this is used to store done_ flags as
        % this module is called once for each session
        sesspath=aas_getsesspath(aap,i,j);
        if (~length(dir(sesspath)))
            [s w]=aas_shell(['mkdir ' sesspath]);
            if (s)
                aas_log(aap,1,sprintf('Problem creating directory for session\n%s',sesspath));
            end
        end
        
        if (iscell(aap.acq_details.subjects(i).seriesnumbers))
            seriesnumbers=aap.acq_details.subjects(i).seriesnumbers{j};
        else
            seriesnumbers=aap.acq_details.subjects(i).seriesnumbers(j);
        end
        
        
        numechoes=length(seriesnumbers);
        if (numechoes>1)
            % Multi-echo
            aas_log(aap,false,sprintf('Session %d has multiple echoes',j));
            echoweightsmoothing=0; % smoothing applied to weights, in mm - 0 to switch off
            % Now taken from XML wrapper
            %             aap.tasklist.currenttask.settings.echoweightmethod='t2star';%'t2star'; % options -  'cnr', 'equal', 't2star'
            subfolder=[];
            fns=[];
            [aap fns DICOMHEADERS]=aas_convertseries_fromstream(aap,i,j,'dicom_epi');
            
            for echoind=1:numechoes
                aap=aas_desc_outputs(aap,i,j,sprintf('epi_echo_%d',echoind),fns{echoind});
            end
            % Now realign first echo to estimate CNR
            allimgs=[];
            for fileind=1:length(fns{1})
                allimgs=strvcat(allimgs,fns{1}{fileind});
            end
            spm_realign(allimgs);
            
            % And copy over the space
            
            minnumimgs=length(fns{1});
            for echoind=2:numechoes
                if (minnumimgs~=length(fns{echoind}))
                    aas_log(aap,false,sprintf('Warning- during multi-echo conversion, not all echoes have the same number of images, using minimum'));
                    minnumimgs=min(minnumimgs,length(fns{echoind}));
                end
                
                for fileind=1:minnumimgs
                    srcimg=fns{1}{fileind};
                    targetimg=fns{echoind}{fileind};
                    M=spm_get_space(srcimg);
                    spm_get_space(targetimg,M);
                end
            end
            
            % Weight according to expected BOLD sensitivity
            Ysnr=[];
            % Get full echo covariance matrix, without loading up all the
            % images at once
            N=length(fns{echoind});
            clear Ytot
            clear Ytotsq
            clear Y1Y2tot
            for fileind=1:N
                for echoind=1:numechoes
                    V=spm_vol(fns{echoind}{fileind});
                    if (fileind==1 & echoind==1)
                        Y={};
                        [Y{1} XYZ]=spm_read_vols(V);
                    else
                        XYZvox=V.mat\[XYZ;ones(1,size(XYZ,2))];
                        Y{echoind}=spm_sample_vol(V,XYZvox(1,:),XYZvox(2,:),XYZvox(3,:),1);
                        Y{echoind}=reshape(Y{echoind},V.dim);
                    end
                end
                for e1=1:numechoes
                    if (fileind==1)
                        Ytot{e1}=Y{e1};
                        Ytotsq{e1}=Y{e1}.^2;
                    else
                        Ytot{e1}=Ytot{e1}+Y{e1};
                        Ytotsq{e1}=Ytotsq{e1}+Y{e1}.^2;
                    end
                    for e2=e1:numechoes
                        if (fileind==1)
                            Y1Y2tot{e1,e2}=Y{e1}.*Y{e2};
                        else
                            Y1Y2tot{e1,e2}=Y1Y2tot{e1,e2}+Y{e1}.*Y{e2};
                        end
                    end
                end
            end
            
            % Covariance CV
            for e1=1:numechoes
                for e2=e1:numechoes
                    CV(e1,e2,:,:,:)=Y1Y2tot{e1,e2}/N-Ytot{e1}.*Ytot{e2}/(N.^2);
                    %./((N*Ytotsq{e1}-Ytot{e1}.^2).*(N*Ytotsq{e2}-Ytot{e2}.^2));
                    CV(e2,e1,:,:,:)=CV(e1,e2,:,:,:);
                end
            end
            
            for e1=1:numechoes
                for e2=e1:numechoes
                    V.dt(1)=spm_type('float32');
                    V.fname=fullfile(sesspath,sprintf('cov_echo_%d_%d.nii',e1,e2));
                    spm_write_vol(V,squeeze(CV(e1,e2,:,:,:)));
                    aap=aas_desc_outputs(aap,i,j,sprintf('cov_echo_%d_%d',e1,e2),V.fname);
                    
                end
            end
            
            % Find temporal SNR
            for echoind=1:numechoes
                Ymean{echoind}=Ytot{echoind}./N;
                Ystd{echoind}=sqrt(N/(N-1)*squeeze(CV(echoind,echoind,:,:,:)));
                Ystd{echoind}(~Ystd{echoind})=nan;
                Ysnr{echoind}=Ymean{echoind}./Ystd{echoind};
                
                %                 dicomheaderfn=fullfile(sesspath,'dicom_headers.mat');
                %                 save(dicomheaderfn,DICOMHEADERS);
                %
                TE{echoind}=DICOMHEADERS{echoind}.EchoTime;
                Ycnr{echoind}=Ysnr{echoind}*TE{echoind};  % Fixed by RC  7/2/2011
                
                
                % Calculate sum for normalisation
                if (echoind==1)
                    Ycnr_tot=Ycnr{1};
                else
                    Ycnr_tot=Ycnr_tot+Ycnr{echoind};
                end
                
                % Write out some files for the record
                subfolder{echoind}=sprintf('echo%d',echoind);
                V.dt(1)=spm_type('float32');
                V.fname=fullfile(sesspath,subfolder{echoind},'echo_mean.nii');
                spm_write_vol(V,Ymean{echoind});
                aap=aas_desc_outputs(aap,i,j,sprintf('epi_echo_%d_mean',echoind),V.fname);
                V.fname=fullfile(sesspath,subfolder{echoind},'echo_std.nii');
                spm_write_vol(V,Ystd{echoind});
                aap=aas_desc_outputs(aap,i,j,sprintf('epi_echo_%d_std',echoind),V.fname);
                V.fname=fullfile(sesspath,subfolder{echoind},'echo_snr.nii');
                spm_write_vol(V,Ysnr{echoind});
                aap=aas_desc_outputs(aap,i,j,sprintf('epi_echo_%d_snr',echoind),V.fname);
                V.fname=fullfile(sesspath,subfolder{echoind},'echo_cnr.nii');
                spm_write_vol(V,Ycnr{echoind});
                aap=aas_desc_outputs(aap,i,j,sprintf('epi_echo_%d_cnr',echoind),V.fname);
            end
            
            for echoind=1:numechoes
                Ymean_SH(echoind,:,:,:) = Ymean{echoind};
                FX(echoind) = TE{echoind};
            end
            FX = FX';
            Ymean_SH(~Ymean_SH)=nan;
            
            FY = log(Ymean_SH);
            C = [FX,ones(size(FY,1),1)]\reshape(FY,size(FY,1),[]);
            C = reshape(C,2,size(FY,2),size(FY,3),size(FY,4));
            rho = squeeze(exp(C(2,:,:,:)));
            T2star = squeeze(1./abs(C(1,:,:,:)));
            % T2star(rho < max(rho(:))/77) = 0; % air
            % T2star(T2star > 2*max(FX))   = 0; % free water
            
            % Write out maps of T2star and proton density
            V.dt(1)=spm_type('float32');
            V.fname=fullfile(sesspath,'t2star.nii');
            spm_write_vol(V,T2star);
            aap=aas_desc_outputs(aap,i,j,'epi_t2star',V.fname);
            
            V.fname=fullfile(sesspath,'protondensity.nii');
            spm_write_vol(V,rho);
            aap=aas_desc_outputs(aap,i,j,'epi_protondensity',V.fname);
            
            % Calculate sum for normalisation
            for echoind=1:numechoes
                if (echoind==1)
                    Ytot=TE{echoind}*exp(-TE{echoind}/T2star);
                else
                    Ytot=Ytot+TE{echoind}*exp(-TE{echoind}/T2star);
                end
            end
            %               % Estimation of covariance hyperparameters
            %             if (strcmp(aap.tasklist.currenttask.settings.echoweightmethod,'covariance'))
            %                 te=[TE{:}];
            %
            %                 t2_min=10;
            %                 t2_max=50;
            %                 t2range=t2_min:t2_max
            %                 numechoessq=numechoes^2;
            %                 numvox=prod(V.dim);
            %                 CV=reshape(CV,[numechoessq numvox]);
            %                 noise=zeros(numechoessq,numvox);
            %                 fit=zeros(3,numvox);
            %                 for t2ind=1:length(t2range)
            %                     t2=t2range(t2ind);
            %                     decay=exp(-1./t2'*te);
            %                     b=repmat(te,[size(decay,1) 1]).*decay;
            %                     X{t2ind}=reshape(b'*b,[numechoessq 1]);
            %                     X{t2ind}(:,2)=reshape(decay'*decay,[numechoessq 1]);
            %                     X{t2ind}(:,3)=reshape(eye(numechoes),[numechoessq 1]);
            %                     XXinv{t2ind}=inv(X{t2ind}'*X{t2ind});
            %                     mask=round(T2star(:))==t2;
            %                     if (any(mask(:)))
            %                         fit=fmincon(@(x) sum(sum((CV(:,mask)-X{t2ind}*x).^2)),ones(3,sum(mask)),-eye(3*sum(mask)),zeros([3*sum(mask) 1]));
            % %                        fit=fmincon(@(x) myfun(x,X{t2ind},CV(:,mask)),ones(3,sum(mask)),-eye(3*sum(mask)),zeros([3*sum(mask) 1]))
            % %lsqlin(X,y,-eye(3),zeros([3 1]))
            %
            % %                        fit(:,mask)=fmincon(X{t2ind},CV(:,mask),-eye(3),zeros([3 1]));
            % %                    fit(:,mask)=XXinv{t2ind}*(X{t2ind}'*CV(:,mask));
            %                     noise(:,mask)=CV(:,mask)-X{t2ind}*fit(:,mask);
            %                     end
            %                 end
            %
            % Measured covariance
            
            
            %                 fit=[covb(:) covi(:) covt(:)]\CV(:,:)
            %                 covfit=fit(1)*covb+fit(2)*covi+fit(3)*covt;
            %               noise=covy-fit(1)*covb;
            
            %             weight=pinv(noise)*b'/sqrt(b*pinv(noise)*b')
            %             weight=weight/sum(weight);
            %
            %             boldest=y*weight;
            %
            
            
            for echoind=numechoes:-1:1
                
                % Choose how to combine different echoes
                switch(aap.tasklist.currenttask.settings.echoweightmethod)
                    case 'cnr'
                        Yweight{echoind}=Ycnr{echoind}./Ycnr_tot;
                    case 'equal'
                        Yweight{echoind}=(1+0.*Ycnr{echoind})/numechoes;
                    case 't2star'
                        Yweight{echoind}=TE{echoind}*exp(-TE{echoind}/T2star)./Ytot;
                    case 'first'
                        Yweight{echoind}=(echoind==1)*ones(size(Ycnr{echoind}));
                    case 'image_and_nulling'
                        Yweight{echoind}=exp(-TE{echoind}/T2star)./Ytot;
                        if (echoind~=1)
                            % Our weight factor for this echo
                            Yweight_nulled{echoind}=TE{echoind}*Yweight{echoind};
                            
                            % Keep track of how much signal each echo will
                            % contribute, by multiplying our weight factor
                            % by the expected image intensity (in Yweight)
                            if (echoind==numechoes)
                                Yweight_tot=Yweight_nulled{echoind}.*Yweight{echoind};
                            else
                                Yweight_tot=Yweight_tot+Yweight_nulled{echoind}.*Yweight{echoind};
                            end
                        else
                            % Weight factor for the first echo, calculated
                            % last (as this loop runs in reverse).
                            % This is designed to cancel out all of the
                            % image signal so far
                            Yweight_nulled{echoind}=-Yweight_tot./Yweight{echoind};
                        end
                end
                
            end
            if (exist('Yweight','var'))
                for echoind=1:numechoes
                    V.fname=fullfile(sesspath,subfolder{echoind},'echo_weight.nii');
                    spm_write_vol(V,Yweight{echoind});
                    if (echoweightsmoothing>0)
                        % smooth the weights
                        [pth nme ext]=fileparts(V.fname);
                        smfn=fullfile(pth,['s' nme ext]);
                        spm_smooth(V.fname,smfn,echoweightsmoothing);
                        V=spm_vol(smfn);
                        Yweight{echoind}=spm_read_vols(V);
                    end


                end
            end
            
            dumpnull=false;
            finalepis={};
            finalepis_nulled={};
            % Now combine the echoes
            te=[TE{:}]';
            % Weighted least squares. Optimise for T2star of 
            % 30ms (grey matter, 3T)
            boldmag=te.*exp(-te./30);
            W=diag(boldmag/sum(boldmag));
            % ridge regression parameter
            k=0.05;
            

            % Now actually do the combination
            % Re-read the images so that we're not outputting the
            % resampled, relaigned ones
            for fileind=1:length(fns{echoind})
                for echoind=1:numechoes
                    V=spm_vol(fns{echoind}{fileind});
                    if (fileind==1) && (echoind==1)
                        Vfirst=V;
                    end
                    Yimg(echoind,:,:,:)=spm_read_vols(V);
                end
                
                if (strcmp(aap.tasklist.currenttask.settings.echoweightmethod,'t2fit'))
                    
                    X=[te./std(te),ones(size(te))];
                    C=inv(X'*W*X +k*eye(2))*X'*W*reshape(log(Yimg),numechoes,[]);
                    C = reshape(C,[2,V.dim]);
                    Ytot = squeeze(exp(C(2,:,:,:)));
                    Ytot_nulled= squeeze(1./abs(C(1,:,:,:)))*std(te);
                    Ytot_nulled(Ytot_nulled>80)=80;
                    Ytot_nulled(Ytot_nulled<5)=5;
                    dumpnull=true;
                else
                    
                    for echoind=1:numechoes
                        Y=squeeze(Yimg(echoind,:,:,:)).*squeeze(Yweight{echoind});
                        
                        switch(aap.tasklist.currenttask.settings.echoweightmethod)
                            case 'image_and_nulling'
                                dumpnull=true;
                                Ynulled=squeeze(Yimg(echoind,:,:,:)).*Yweight_nulled{echoind};
                                if (echoind==1)
                                    Ytot_nulled=Ynulled;
                                else
                                    Ytot_nulled=Ytot_nulled+Ynulled;
                                end
                                
                        end
                        if (echoind==1)
                            Ytot=Y;
                        else
                            Ytot=Ytot+Y;
                        end
                    end
                end
                
                % Output combined file
                [pth fle ext]=fileparts(fns{echoind}{fileind});
                Vfirst.fname=fullfile(sesspath,[fle '_multiechocombo' ext]);
                finalepis{fileind}=Vfirst.fname;
                spm_write_vol(Vfirst,Ytot);
                
                if (dumpnull)
                    [pth fle ext]=fileparts(fns{echoind}{fileind});
                    Vfirst.fname=fullfile(sesspath,[fle '_multiechocombo_nulled' ext]);
                    finalepis_nulled{fileind}=Vfirst.fname;
                    spm_write_vol(Vfirst,Ytot_nulled);
                end
            end
            
        else
            
            %% Single echo code
            subfolder=[];
            fns=[];
            
            % No output stream, so won't save (it is saved below with dummy
            % separation)
            [aap fns DICOMHEADERS]=aas_convertseries_fromstream(aap,i,j,'dicom_epi');
            finalepis=fns;
            Ysnr=[];
            % Find temporal SNR
            for fileind=1:length(fns)
                V=spm_vol(fns{fileind});
                if (fileind==1)
                    [Y XYZ]=spm_read_vols(V);
                    Ytot=Y;
                    Ytotsq=Y.*Y;
                else
                    XYZvox=V.mat\[XYZ;ones(1,size(XYZ,2))];
                    Y=spm_sample_vol(V,XYZvox(1,:),XYZvox(2,:),XYZvox(3,:),1);
                    Y=reshape(Y,size(Ytot));
                    Ytot=Ytot+Y;
                    Ytotsq=Ytotsq+Y.*Y;
                end
            end
            N=length(fns);
            Ymean=Ytot./N;
            Ystd=sqrt(N/(N-1)*(Ytotsq/N-Ymean.^2));
            % Sometimes when the value is close to zero, we get slight
            % negative values, which render some voxels "imaginary"
            Ystd(imag(Ystd)~=0) = 0; % Remove these [AVG]
            Ystd(~Ystd)=nan;
            Ysnr=Ymean./Ystd;
            
            TE=DICOMHEADERS{1}.EchoTime;
            Ycnr=Ysnr*TE;
            
            % Write out some files for the record
            V.dt(1)=spm_type('float32');
            V.fname=fullfile(sesspath,'echo_mean.nii');
            spm_write_vol(V,Ymean);
            aap=aas_desc_outputs(aap,i,j,'epi_mean',V.fname);
            V.fname=fullfile(sesspath,'echo_std.nii');
            spm_write_vol(V,Ystd);
            aap=aas_desc_outputs(aap,i,j,'epi_std',V.fname);
            V.fname=fullfile(sesspath,'echo_snr.nii');
            spm_write_vol(V,Ysnr);
            aap=aas_desc_outputs(aap,i,j,'epi_snr',V.fname);
            V.fname=fullfile(sesspath,'echo_cnr.nii');
            spm_write_vol(V,Ycnr);
            aap=aas_desc_outputs(aap,i,j,'epi_cnr',V.fname);
            
            boldmap=TE*Ymean;
            V.fname=fullfile(sesspath,'boldmap.nii');
            spm_write_vol(V,boldmap);
            aap=aas_desc_outputs(aap,i,j,'epi_boldmap',V.fname);
            
            dumpnull=false;
        end
        
        % Now move dummy scans to dummy_scans directory
        dummypath=fullfile(sesspath,'dummy_scans');
        aap=aas_makedir(aap,dummypath);
        
        dummylist=[];
        for d=1:aap.acq_details.numdummies
            cmd=['mv ' finalepis{d} ' ' dummypath];
            [pth nme ext]=fileparts(finalepis{d});
            dummylist=strvcat(dummylist,fullfile('dummy_scans',[nme ext]));
            [s w]=aas_shell(cmd);
            if (s)
                aas_log(aap,1,sprintf('Problem moving dummy scan\n%s\nto\n%s\n',convertedfns{d},dummypath));
            end
        end
        
        % And describe outputs
        aap=aas_desc_outputs(aap,i,j,'dummyscans',dummylist);
        dcmhdrfn=fullfile(sesspath,'dicom_headers.mat');
        save(dcmhdrfn,'DICOMHEADERS');
        aap=aas_desc_outputs(aap,i,j,'epi',finalepis(aap.acq_details.numdummies+1:end));
        
        if dumpnull
            aap=aas_desc_outputs(aap,i,j,'dummy_scans_nulled',finalepis_nulled(1:aap.acq_details.numdummies));
            aap=aas_desc_outputs(aap,i,j,'epi_nulled',finalepis_nulled(aap.acq_details.numdummies+1:end));
        end
        
        aap=aas_desc_outputs(aap,i,j,'epi_header',dcmhdrfn);
        
        
        % And describe outputs
        aap=aas_desc_outputs(aap,i,j,'dummyscans',dummylist);
        dcmhdrfn=fullfile(sesspath,'dicom_headers.mat');
        save(dcmhdrfn,'DICOMHEADERS');
        aap=aas_desc_outputs(aap,i,j,'epi_header',dcmhdrfn);
        
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end
