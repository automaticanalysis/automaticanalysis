% AA module - first level statistics for PPI
% [aap,resp]=aamod_ppi_model(aap,task,i)
% You will need to make a local copy of this module and then modify it to
%  reflect the particular experiment you conducted
% First-level model based on original by Adam Hampshire MRC CBU Cambridge Feb 2006
% Modified for aa by Rhodri Cusack MRC CBU Mar-May 2006

function [aap,resp]=aamod_ppi_model(aap,task,i)

resp='';

switch task
    case 'domain'
        resp='subject';   % this module needs to be run once per subject

    case 'description'
        resp='SPM5 PPI model';

    case 'summary'
        subjpath=aas_getsubjpath(i);
        resp=sprintf('SPM5 PPI model for %s\n',subjpath);

    case 'report'

    case 'doit'

        % YOU WILL NEED TO CHANGE THIS SECTION
        % condition names
        cnames={'small_b_arm_large_b_arm','small_b_arm_large_nb_arm','small_b_leg_large_b_leg','small_b_leg_large_nb_leg','small_move_arm_left_large_move_arm_left','small_move_arm_right_large_move_arm_right','small_move_leg_left_large_move_leg_left','small_move_leg_right_large_move_leg_right','small_nb_arm_large_b_arm','small_nb_arm_large_nb_arm','small_nb_leg_large_b_leg','small_nb_leg_large_nb_leg','response'};
        ncond=length(cnames);
        % high pass filter
        hpf = 128;


        % include movement parameters?
        incmoves = 1;


        %get subject directory
        cwd=pwd;
        subj_dir = aas_getsubjpath(aap,i);

        %number of blocks
        nsess = length(aap.acq_details.sessions);


        % retrieve TR from DICOM header
        if (length(aap.spmanalysis.TRs)==0)
            DICOMHEADERS=load(fullfile(aas_getsesspath(aap,i,1),'dicom_headers'));
            aap.spmanalysis.TRs=DICOMHEADERS.DICOMHEADERS{1}.RepetitionTime/1000;
        end;



        modeldur = 0;

        movefilt = '^rp_.*\.txt$';
        mnames = {'x_trans' 'y_trans' 'z_trans' 'x_rot' 'y_rot' 'z_rot'};



        %-----------------------------------------------------------------------
        %Design setup
        %-----------------------------------------------------------------------
        clear SPM

        SPM.xY.RT = aap.spmanalysis.TRs;
        SPM.xGX.iGXcalc = 'None';
        SPM.xVi.form = 'none'; %AR(1)';

        % not used but set up values anyway or SPM prompts!
        SPM.xBF.T          = 16;                % number of time bins per scan
        SPM.xBF.T0         = 1;                 % first time bin (see slice timing) - middle of TA
        SPM.xBF.UNITS      = 'scans';           % OPTIONS: 'scans'|'secs' for onsets
        SPM.xBF.Volterra   = 1;                 % OPTIONS: 1|2 = order of convolution
        SPM.xBF.name       = 'hrf';
        SPM.xBF.length     = 32;              % length in seconds
        SPM.xBF.order      = 1;                 % order of basis set

        % Pre- calculate some things that are needed for getting the model
        % columns for the PPI
        SPM.xBF.dt = SPM.xY.RT/SPM.xBF.T;   % calc dt
        % Get basis functions
        try
            SPM.xBF.bf;
        catch
            SPM.xBF = spm_get_bf(SPM.xBF);
        end

        subdata = aas_getsubjpath(aap,i);
        % New option to allow suffix to output file in extraparameters
        if (isfield(aap.tasklist.currenttask.extraparameters,'stats_suffix'))
            stats_suffix=aap.tasklist.currenttask.extraparameters.stats_suffix;
        else
            stats_suffix=[];
        end;

        anadir = fullfile(subdata,[aap.directory_conventions.stats_singlesubj stats_suffix]);
        if exist(anadir)~=7; mkdir(subdata,[aap.directory_conventions.stats_singlesubj stats_suffix]);end
        cd(anadir);

        % load up VOI data
        voidir=fullfile(aas_getsubjpath(aap,i),aap.directory_conventions.voidirname);
        voifn=fullfile(voidir,aap.directory_conventions.voifilename);
        vois=load(voifn);
        vois=vois.vois;
        voinames=[];
        for r=1:length(vois)
            voinames=[voinames {vois{r}.name}];
        end;

        tc = 0;
        allfiles='';
        for sess = 1:nsess
            tc = tc+1;
            files = aas_getimages(aap,i,sess,aap.tasklist.currenttask.epiprefix);
            evdir = aas_getsesspath(aap,i,sess);

            cind=1;
            for c = 1:ncond
                % YOU MAY NEED TO CHANGE THIS
                % It currently expects filenames
                % ons_[sessionnumber]_[conditionname].txt (e.g., ons_1_flash.txt) for event onsets (in scans)
                % and dur_[sessionnumber]_[conditionname].txt for event durations (in scans)

                efile = fullfile(evdir,['ons_' cnames{c} '.txt']);
                if (exist(efile,'file'))
                    ons = load(efile)
                    efile2 = fullfile(evdir,['dur_'  cnames{c} '.txt']);
                    durs = load(efile2)

                    for isons=0:0 %1
                        if (isons==0)
                            colname=cnames(c);
                        else
                            colname=[cnames(c) '_onset'];
                        end;
                        SPM.Sess(tc).U(cind) = struct(...
                            'ons',ons,...
                            'dur',durs,...
                            'name',{colname},...
                            'P',struct('name','none'))
                        cind=cind+1;
                        durs=zeros(size(durs));
                    end;
                end;
            end
            SPM.Sess(tc).C.C = [];
            SPM.Sess(tc).C.name = [];
            SPM.nscan(tc) = size(files,1);

            SPM.xX.K(tc).HParam = hpf;

            allfiles = strvcat(allfiles,files);

        end;

        cons=[];
        cd (anadir)

        SPM.xY.P = allfiles;




        % Now, for each ppi required
        for p=1:length(aap.spmanalysis.ppis)
            % Get VOI data...
            voinum=find(strcmp(voinames,aap.spmanalysis.ppis{p}.voiname));
            if (isempty(voinum))
                aas_log(aap,true,sprintf('Cannot find VOI %s as specified in PPI %s\n',aap.spmanalysis.ppis{p}.voiname,aap.spmanalysis.ppis{p}.ppiname));
            end;

            % Now generate PPI column for each session...
            for sess=1:nsess
                % PPI PART 1: data in VOI
                % Filter VOI data (no whitening at present)
                K=[];
                K(1).row=[1:SPM.nscan(sess)];
                K(1).HParam=hpf;
                K(1).RT=aap.spmanalysis.TRs;
                y=vois{voinum}.data{sess}.raw;
                y=spm_filter(K,y);
                % ...compute regional response in terms of first eigenvariate
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

                % PPI PART 2: make context regressor
                % Create unfiltered design matrix
                U=spm_get_ons(SPM,sess);
                X=spm_Volterra(U,SPM.xBF.bf,1);
                X = X([0:(SPM.nscan(sess) - 1)]*SPM.xBF.T + SPM.xBF.T0 + 32,:);
                % ...make context variable
                context=X*(aap.spmanalysis.ppis{p}.contrast(:));
                % ...multiply VOI data for this session by context variable

                % PPI COMBINE THE TWO PARTS
                Y=Y.*context;
                % ...and add to regressors
                SPM.Sess(sess).C.C    = [SPM.Sess(sess).C.C Y];
                SPM.Sess(sess).C.name = [SPM.Sess(sess).C.name {aap.spmanalysis.ppis{p}.ppiname}]; % [1 x c cell]   names
            end;
        end;

        % And add movement parameters on the end
        for sess=1:nsess
            if incmoves==1
                sessdata=aas_getsesspath(aap,i,sess);
                mfname = spm_select ('List', sessdata, movefilt);
                moves = load(fullfile(sessdata,mfname));
                SPM.Sess(sess).C.C    = [SPM.Sess(sess).C.C moves];     % [n x c double] covariates
                SPM.Sess(sess).C.name = [SPM.Sess(sess).C.name mnames]; % [1 x c cell]   names
            end
        end;

        % Set up design
        SPMdes = spm_fmri_spm_ui(SPM);

        % This turns off masking around the brain
        %        SPMdes.xM.TH=-inf(size(SPMdes.xM.TH));

        spm_unlink(fullfile('.', 'mask.img')); % avoid overwrite dialog

        % now check real covariates and nuisance variables are
        % specified correctly
        for sess=1:nsess
            sesscol=(sess-1)*7+1;
            SPMdes.xX.iG=[SPMdes.xX.iG (sesscol+1):(sesscol+6)];
        end;
        SPMdes.xX.iC=setdiff(SPMdes.xX.iC,SPMdes.xX.iG);
        SPM=SPMdes;
        save SPM SPM
        % and estimate..
        SPMest = spm_spm(SPMdes);

        cd (cwd);


    case 'checkrequirements'

    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end;



