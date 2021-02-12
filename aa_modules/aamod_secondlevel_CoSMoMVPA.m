function [aap,resp]=aamod_secondlevel_CoSMoMVPA(aap,task,subj)
resp='';

switch task
    case 'report'
        f{1} = fullfile(aas_getstudypath(aap),'diagnostic_overlay_3_001.jpg');
        if exist(f{1},'file')
            f{2} = fullfile(aas_getstudypath(aap),'diagnostic_render.jpg');

            tstat = dlmread(strrep(f{1},'_overlay_3_001.jpg','.txt'));
            aap = aas_report_add(aap,[],sprintf('T = %2.2f - %2.2f</tr><tr>',tstat(1),tstat(2)));
            
            aap = aas_report_add(aap,[],'<table><tr>');
            for i = 1:2
                aap = aas_report_add(aap,[],'<td>');
                aap=aas_report_addimage(aap,[],f{i});
                aap = aas_report_add(aap,[],'</td>');
            end
            aap = aas_report_add(aap,[],'</tr></table>');
        end
    case 'doit'
        %% Prepare
        % Initialise Cosmo
        [junk, MVPA] = aas_cache_get(aap,'cosmomvpa');
        MVPA.load;
        
        cosmodir = fullfile(aas_getstudypath(aap),'cosmo');
        aas_makedir(aap,cosmodir);
        
        % Background
        switch  aap.tasklist.currenttask.settings.overlay.template
            case 'averaged_structurals'
                tmpfile = [];
                inpstreams = aas_getstreams(aap,'input');
                if aas_stream_has_contents(aap,inpstreams{end})
                    for s = 1:numel(aap.acq_details.subjects)
                        tmpfile{s} = aas_getfiles_bystream(aap,s,inpstreams{end});
                    end
                    aas_makedir(aap,fullfile(aas_getstudypath(aap),'structural'));
                    tmpfile = spm_imcalc(char(tmpfile),fullfile(aas_getstudypath(aap),'structural','background.nii'),'mean(X)',struct('dmtx',1,'mask',1));
                    tmpfile = tmpfile.fname;
                else
                    aas_log(aap, true, sprintf('%s: No structural stream found. Define one, or use SPMT1 template option. Exiting...', mfilename));
                end
                
            case 'SPMT1'
                % assume a reasonable default location, but assume the user put
                % the correct location in aap.dir_con.T1template if it's not empty
                tmpfile = 'toolbox/OldNorm/T1.nii';
                if ~isempty(aap.directory_conventions.T1template) tmpfile = aap.directory_conventions.T1template; end
                if (tmpfile(1) ~= '/'), tmpfile = fullfile(fileparts(which('spm')),tmpfile); end
                if ~exist(tmpfile,'file')
                    aas_log(aap, true, sprintf('%s: SPM T1 template not found. Exiting...', mfilename));
                end
                
            otherwise
                aas_log(aap, true, sprintf('%s: Unknown template option. Exiting...', mfilename));
        end
        
        % maps
        mapfiles = {};
        for subj = 1:numel(aap.acq_details.subjects)
            mapfiles{subj} = aas_getfiles_bystream(aap,'subject',subj,aas_getstreams(aap,'input',1));
        end
        spm_file_merge(mapfiles,fullfile(cosmodir,'group.nii'));
        
        tol = 10e-3;
        V = spm_vol(fullfile(cosmodir,'group.nii'));
        Y = spm_read_vols(V);
        M = Y>tol;
        M = min(M,[],4);
        V = V(1);
        V.fname = spm_file(V.fname,'basename','mask');
        V.descrip = 'mask';
        spm_write_vol(V,M);

        params.niter = aas_getsetting(aap,'inference.iteration');
        
        % model
        % - one-sample
        groupmodel = ones(numel(mapfiles),1);
        subjectmodel = (1:numel(mapfiles))';
        switch aas_getstreams(aap,'input',1)
            case 'RSAmap'
                params.h0_mean = 0;
            case 'Cmap'
                k = aas_getsourcestage(aap,'aamod_CoSMoMVPA','Cmap');
                srcaap = aas_setcurrenttask(aap,k);
                params.h0_mean = 1/numel(aas_getsetting(srcaap,'itemList'));
        end
        
        if aas_getsetting(aap,'numworker') > 1, params.nproc = aas_getsetting(aap,'numworker'); end
        
        %% Run
        g_ds=cosmo_fmri_dataset(fullfile(cosmodir,'group.nii'),'mask',fullfile(cosmodir,'mask.nii'),'targets',groupmodel,'chunks',subjectmodel);
        cluster_nbrhood=cosmo_cluster_neighborhood(g_ds); 
        
        params.cluster_stat = aas_getsetting(aap,'inference.correction');
        switch aas_getsetting(aap,'inference.method')
            case 'montecarlo'
                if ~strcmp(aas_getsetting(aap,'inference.correction'),'tfce')
                    params.p_uncorrected = aas_getsetting(aap,'inference.pclusterforming');
                end
        end
        stat_map=cosmo_montecarlo_cluster_stat(g_ds, cluster_nbrhood, params);
                    
        %% Postprocess
        if aas_getsetting(aap,'numworker') > 1, delete(gcp('nocreate')); end
        
        cosmo_map2fmri(stat_map,fullfile(cosmodir,'zmap.nii'));

        p_map = stat_map;
        p_map.samples = 2*normcdf(stat_map.samples);
        cosmo_map2fmri(p_map,fullfile(cosmodir,'pmap.nii'));
        sel = p_map.samples < aas_getsetting(aap,'inference.p');
        p_map.samples = p_map.samples(sel);
        p_map.fa.i = p_map.fa.i(sel);
        p_map.fa.j = p_map.fa.j(sel);
        p_map.fa.k = p_map.fa.k(sel);
        cosmo_map2fmri(p_map,fullfile(cosmodir,'thresh_pmap.nii'));
        %map_overlay(tmpfile,{fullfile(cosmodir,'thresh_pmap.nii')},'axial',-60:6:60);

        %% Cleanup
        MVPA.unload;
            
    case 'checkrequirements'
        if ~aas_cache_get(aap,'cosmomvpa'), aas_log(aap,true,'CoSMoMVPA is not found'); end
       
end
end