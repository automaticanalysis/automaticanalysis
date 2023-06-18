function [aap,resp] = aamod_render(aap, task)
%
% generate jpeg renders
%
% Change History
%
% spring 2023 [MSJ] -- new
%

resp='';

switch task
	
    case 'report'
        
    case 'doit'	
       
        
      % 1) setup ==========================================================
  
      switch (aap.tasklist.currenttask.settings.renderer)
          
      case 'spm'
          
        % reference:
        %     
        % overlay_nifti(nii_fname, template_fname, render_fname, savefile_fname, image_label, brightness)
      
        if ~isempty(aap.tasklist.currenttask.settings.spm_render.template_fname)
            template_fname = aap.tasklist.currenttask.settings.spm_render.template_fname;
        else
            template_fname = aap.directory_conventions.T1template;
        end
        
        if ~isempty(aap.tasklist.currenttask.settings.spm_render.render_fname)
            render_fname = aap.tasklist.currenttask.settings.spm_render.render_fname;
        else
            render_fname = aap.directory_conventions.Render;
        end
        
        image_label = aap.tasklist.currenttask.settings.spm_render.label;
        brightness = aap.tasklist.currenttask.settings.spm_render.brightness;
        
         
      case 'neurodot'
                
        [ ~, NEURODOT ] = aas_cache_get(aap,'neurodot');
        NEURODOT.load;

        % Loading underlays
        
        % these are in the NeuroDOT folder:
        % NeuroDOT_NITRC-Release/NeuroDOT/Support_Files/

        load('MNI_coord_meshes_164k', 'MNIl', 'MNIr');

        [ ~,infoMNI] = LoadVolumetricData('mni152nl_T1_on_333_nifti', '', 'nii');

        params = [];
        
        params.scale   = aap.tasklist.currenttask.settings.neurodot_render.scale;
        params.pthresh = aap.tasklist.currenttask.settings.neurodot_render.positive_threshold;
        params.nthresh = aap.tasklist.currenttask.settings.neurodot_render.negative_threshold;

        params.Cmap = aap.tasklist.currenttask.settings.neurodot_render.colormap;
        params.view = aap.tasklist.currenttask.settings.neurodot_render.view;
        params.CBar_on = aap.tasklist.currenttask.settings.neurodot_render.show_colorbar;
                 
      case 'washu'

        cfg = [];
        
        % predefine some useful layouts with reasonable defaults

        cfg.inflate = aap.tasklist.currenttask.settings.washu_render.inflate;

        switch (aap.tasklist.currenttask.settings.washu_render.layout)

            case '4view'

                cfg.plots = [1 2 3 4]; % L,R,T,B
                cfg.plot1pos = [.05 .4 .35 .3];
                cfg.plot2pos = [.6 .4 .35 .3];
                cfg.plot4pos = [.351 .18 .3 .35];   % T on top
                cfg.plot3pos = [.350 .55 .3 .35];  

            case 'lefthemi'

                cfg.plots = [1]; % L
                cfg.plot1pos = [.1 .1 .8 .8];

            case 'righthemi'

                cfg.plots = [2]; % R
                cfg.plot1pos = [.1 .1 .8 .8];

            otherwise

                aas_log(aap, true, sprintf('Unrecognized washu layout option %s. Halting...', options));

        end

        % override defaults with any custom settings provided

        if (~isempty(aap.tasklist.currenttask.settings.washu_render.cfg))

            cfg_from_xml = aap.tasklist.currenttask.settings.washu_render.cfg;

            if (isfield(cfg_from_xml,'plots'))
                cfg.plots = cfg_from_xml.plots;
            end

            if (isfield(cfg_from_xml,'plot1pos'))
                cfg.plot1pos = cfg_from_xml.plot1pos;
            end

            if (isfield(cfg_from_xml,'plot2pos'))
                cfg.plot2pos = cfg_from_xml.plot2pos;
            end

            if (isfield(cfg_from_xml,'plot3pos'))
                cfg.plot3pos = cfg_from_xml.plot3pos;
            end

            if (isfield(cfg_from_xml,'plot4pos'))
                cfg.plot4pos = cfg_from_xml.plot4pos;
            end

            if (isfield(cfg_from_xml,'inflate'))
                cfg.inflate = cfg_from_xml.inflate;
            end
            
            if (isfield(cfg_from_xml,'colorscale'))
                cfg.colorscale = cfg_from_xml.colorscale;
            end
          

        end
                 
        otherwise
          
          aas_log(aap,'true',sprintf('Unknown renderer %s. Halting.', aap.tasklist.currenttask.settings.renderer));      
        
        end % switch renderer for setup
        
        % stream input ==========================================
       
        map_inputstream_struct = aap.tasklist.currenttask.inputstreams(1).stream{1};

        if (~aas_stream_has_contents(aap, map_inputstream_struct.CONTENT))
            aas_log(aap, true, sprintf('No data found. Halting...'));
        end

        map_fname_list = cellstr(aas_getfiles_bystream(aap, map_inputstream_struct.CONTENT));
            
        % if we are to apply a common colorscale across all renders,
        % we need to do an initial pass over the data and compute
        % global min/max values
       
        if (aap.tasklist.currenttask.settings.use_common_colorscale)

            all_min = inf;
            all_max = -inf;

            for findex = 1:numel(map_fname_list)

                header = spm_vol(map_fname_list{findex});
                data = spm_read_vols(header);
                all_min = min([all_min;data(:)],[],'omitnan');
                all_max = max([all_max;data(:)],[],'omitnan');

            end
                        
            % we use these global max/min values to scale the render below
            % (instead of computing this_min and this_max for each file)
            
            this_min = all_min;
            this_max = all_max;
            
            recompute_colorscale = false;
            
        else
            
            recompute_colorscale = true;

        end      
 
        
        % 2) render ==============================================================         

        % there are potentially many renderfiles -- put in a "render" directory

        savedir_name = fullfile(aas_getstudypath(aap),'renders');

        if (~exist(savedir_name,'dir') && ~mkdir(savedir_name))
            aas_log(aap, true, sprintf('Could not create save directory %s. Halting...',savedir_name));
        end
        
        % loop over inputs and render/save each image

        for findex = 1:numel(map_fname_list)           
                       
            if (recompute_colorscale)
            
                header = spm_vol(map_fname_list{findex});
                data = spm_read_vols(header);

                this_min = min(data(:),[],'omitnan');
                this_max = max(data(:),[],'omitnan');
                 
            end
             
            nii_fname = map_fname_list{findex};
            [~,n,~] = fileparts(nii_fname);

            jpg_fname = fullfile(savedir_name,[n '.jpg']);


            switch aap.tasklist.currenttask.settings.renderer

                case 'spm'

                    % overlay_nifti is in aa/extrafunctions...

                    % template_fname, render_fname, image_label, brightness set in part (1)
                                        
                    [ errflag,errstring ] = overlay_nifti(nii_fname, template_fname, render_fname, jpg_fname, image_label, brightness);

                    if (errflag)
                        aas_log(aap,'true',sprintf('SPM render returned error %s. Halting.', errstring));          
                    end
                    
                case 'neurodot'

                    % reference
                    % nii_overlay = LoadVolumetricData('SUB01SESS01_SUB02SESS01','','nii');
                    
                    [p,n,~] = fileparts(map_fname_list{findex});
                    nii_overlay = LoadVolumetricData(fullfile(p,n),'','nii');
                    
                    % note "this_max" was set above (either for this file)
                    % or earlier as the global max over all files)
                    
                    params.Scale =  params.scale * this_max;
                    params.Th.P = params.pthresh * this_max;
                    params.Th.N = params.nthresh * this_max;
                 
                    if (aap.tasklist.currenttask.settings.plot_positive_only)
                        params.Th.N = 0;
                    end

                    if (aap.tasklist.currenttask.settings.plot_negative_only)
                        params.Th.P = 0;
                    end

                    PlotInterpSurfMesh(nii_overlay, MNIl, MNIr, infoMNI, params);

                    saveas(gcf,jpg_fname,'jpeg')
                    close(gcf);


                case 'washu'
                   
                    % explicitly defining a colorscale in the cfg struct
                    % overrides all other colorscale options

                    if (~isfield(cfg,'colorscale'))
                        
                        % explicitly defining a colorscale in the xml
                        % overrides dynamic scaling

                        if ~isempty(aap.tasklist.currenttask.settings.washu_render.colorscale)

                            cfg.colorscale = aap.tasklist.currenttask.settings.washu_render.colorscale;

                        else

                            % note "this_min" and "this_max" were set above (either 
                            % for this file or as the global values over all files)

                            cfg.colorscale = [ this_min this_max ];

                            if (aap.tasklist.currenttask.settings.plot_positive_only)
                                cfg.colorscale = [ 0 this_max ]; 
                            end

                            if (aap.tasklist.currenttask.settings.plot_negative_only)
                                cfg.colorscale = [ this_min 0 ];
                            end


                        end
                        
                    end
                    
                    washu_surfacerender(map_fname_list{findex}, [], cfg, jpg_fname);

            end % switch renderer

        end % loop over files

        
        % 3) cleanup ======================================================
        
        switch aap.tasklist.currenttask.settings.renderer
               
            case 'spm'
                
            case 'washu'
                
            case 'neurodot'

                NEURODOT.unload;
        end
        
        
    case 'checkrequirements'
        
        if (strcmp(aap.tasklist.currenttask.settings.renderer,'neurodot'))

            if ~aas_cache_get(aap,'neurodot')
                aas_log(aap,true,'NeuroDot Toolbox not found (required to use neurodot renderer).'); 
            end

        end
        
    otherwise
        aas_log(aap,true,sprintf('Unknown task %s', task));
        

end % task switch

end









