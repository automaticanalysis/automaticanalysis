% AA module - movie correlations
% i=subject num; j=session
% Rhodri Cusack MRC CBU 2010

function [aap,resp]=aamod_moviecorr_summary(aap,task,isc_sess)

resp='';

switch task
    case 'report'
        
    case 'doit'
        
        statspth=aas_getpath_bydomain(aap,'isc_session',isc_sess);
        nsub=length(aap.acq_details.subjects);

        clear SPM
        
        %-Assemble SPM structure
        %=======================================================================
        
        SPM.nscan = nsub;

        for i=1:nsub
            
            img=aas_getfiles_bystream(aap,'isc_subject',[isc_sess i],'moviecorr_loo');
            V=spm_vol(img);
            Y=spm_read_vols(V);
            if i==1 
                Ytot=Y;
            else
                Ytot=Ytot+Y;
            end;
            
            SPM.xY.P{i}   = img;
            SPM.xY.VY(i)   = spm_vol(SPM.xY.P{i});
        end;
        
        
        % Write mean moviecorr and describe as output
        Ytot=Ytot/nsub;
        V.fname=fullfile(statspth,['mean_moviecorr_loo.nii']);
        spm_write_vol(V,Ytot);
        aap=aas_desc_outputs(aap,'isc_session',isc_sess,'mean_moviecorr_loo',V.fname);
        
        
        % Build SPM model
        
        SPM.xX = struct(	'X',	ones(nsub,1),...
            'iH',1,'iC',zeros(1,0),'iB',zeros(1,0),'iG',zeros(1,0),...
            'name',{{'mean'}},'I',[[1:nsub]' ones(nsub,3)],...
            'sF',{{'obs'  ''  ''  ''}});
        
        SPM.xC = [];
        
        SPM.xGX = struct(...
            'iGXcalc',1,	'sGXcalc','omit',				'rg',[],...
            'iGMsca',9,	'sGMsca','<no grand Mean scaling>',...
            'GM',0,		'gSF',ones(nsub,1),...
            'iGC',	12,	'sGC',	'(redundant: not doing AnCova)',	'gc',[],...
            'iGloNorm',9,	'sGloNorm','<no global normalisation>');
        
        SPM.xVi	= struct('iid',1,'V',speye(nsub));
        
        Mdes 	= struct(	'Analysis_threshold',	{'None (-Inf)'},...
            'Implicit_masking',	{'Yes: NaNs treated as missing'},...
            'Explicit_masking',	{'Yes: SPM2 Brain Mask'});
        
        %SPM.xM	= struct('T',-Inf,'TH',ones(nsub*2,1)*-Inf,...
        %		 'I',1,'VM',spm_vol('/home/rh01/SPM/spm5/apriori/brainmask.nii'),'xs',Mdes);
        
        SPM.xM	= struct('T',-Inf,'TH',ones(nsub*2,1)*-Inf,...
            'I',1,'VM',[],'xs',Mdes);
        
        Pdes 	= {{'1 condition, +0 covariate, +0 block, +0 nuisance'; '1 total, having 1 degrees of freedom'; 'leaving 8 degrees of freedom from 9 images'}};
        
        SPM.xsDes = struct(	'Design',		{'One sample t-test'},...
            'Global_calculation',	{'omit'},...
            'Grand_mean_scaling',	{'<no grand Mean scaling>'},...
            'Global_normalisation',	{'<no global normalisation>'},...
            'Parameters',		Pdes);
        
        % Estimate parameters
        %===========================================================================
        cwd=pwd;
        cd(statspth);
        spm_unlink(fullfile('.', 'mask.img')); % avoid overwrite dialog
        rSPM = spm_spm(SPM);
        cd(cwd)
        
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end;














