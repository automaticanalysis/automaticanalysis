% AA module - movie correlations
% i=subject num; j=session
% Rhodri Cusack MRC CBU 2010

function [aap,resp]=aamod_moviecorr_summary(aap,task)

resp='';

switch task
    case 'report'
        
    case 'doit'
        conlist={{'session','moviecorr_meantimecourse'}};
        cwd=pwd;
        for c=1:length(conlist)
            nsub=length(aap.acq_details.subjects);
            nsess=length(aap.acq_details.selected_sessions);
            for j=aap.acq_details.selected_sessions
                condir=fullfile(aas_getstudypath(aap),sprintf('%s_%d',conlist{c}{2},j));
                aap=aas_makedir(aap,condir);
                cd(condir);
                clear SPM
                
                %-Assemble SPM structure
                %=======================================================================
                
                SPM.nscan = nsub;
                for i=1:nsub
                    switch (conlist{c}{1})
                        case 'subject'
                            img=aas_getfiles_bystream(aap,i,conlist{c}{2});
                        case 'session'
                            img=aas_getfiles_bystream(aap,i,j,conlist{c}{2});
                    end;
                    V=spm_vol(img);
                    Y=spm_read_vols(V);
                    if (i==1 && j==aap.acq_details.selected_sessions(1))
                        Ytot=Y;
                    else
                        Ytot=Ytot+Y;
                    end;
                    
                    SPM.xY.P{i}   = img;
                    SPM.xY.VY(i)   = spm_vol(SPM.xY.P{i});
                end;
                
                
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
                spm_unlink(fullfile('.', 'mask.img')); % avoid overwrite dialog
                rSPM = spm_spm(SPM);
                
            end;
            
            Ytot=Ytot./(nsub*nsess);
            V.fname=fullfile(aas_getstudypath(aap),['mean' conlist{c}{2} '.nii']);
            spm_write_vol(V,Ytot);
        end;
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end;














