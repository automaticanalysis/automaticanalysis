function [aap,resp] = aamod_lossplot(aap, task)
%
% generate a "lossplot" (threshold v. %dataloss across all sessions / all
% subjects) for a given metric (supplied in the input stream metric_data)
%
% the expected format of metric_data goes metric_data.name, for example
% metric_data.FD, metric_data.DVARS, etc.
%
% lossplots for all metrics defined in metric_data are saved. Also,
% values of the metric for TOI (thresholds of interest) are saved in
% the output stream metric_thresholds. The format mimics the format
% of metric_data (e.g., metric_data.FD => metric_thresholds.FD.onepercent,
% metric_thresholds.FD.twoprecent, etc).
%

resp='';

switch task
	
    case 'report'
		
        fdiag = dir(fullfile(aap.acq_details.root,'*.jpg'));
        for d = 1:numel(fdiag)
            aap = aas_report_add(aap,[],'<table><tr><td>');
            aap=aas_report_addimage(aap,[],fullfile(aap.acq_details.root,fdiag(d).name));
            aap = aas_report_add(aap,[],'</td></tr></table>');
		end 
		
    case 'doit'
														
		% we get metric_names as the fieldnames of metric_data,
		% which specifies the the time series data to use
		
		% this is a bit clunky, but we need the number of metrics we're going to be plotting
		% in order to loop over them - try subject 1, first session (remember - we're a study level module):
		
		temp = load(aas_getfiles_bystream(aap, 1, min(aap.acq_details.selected_sessions), 'metric_data'));
		metric_names = fieldnames(temp.metric_data);
	
		TOI = [ 1.0 2.0 5.0 10.0 20.0 40.0 ];	% dataloss thresholds of interest
		quantiles_to_plot = [0:0.01:1];			% this needs good coverage of transition zone (ergo 0.01)
		
		nsub = length(aap.acq_details.subjects);
			
		for mindex = 1:numel(metric_names)
			
			metric_name = metric_names{mindex};
			
			metric_subject_values = zeros(nsub,length(quantiles_to_plot));
			
			% gather data from each subject across all sessions for this metric
			% this is inefficient (in that we reread the metric files on each
			% pass) but it makes for simpler code.
			
			% the movegui trick here centers the window, which looks cool
			% but it breaks plotting if running headless (i.e. cluster)

			if (strcmp(aap.options.wheretoprocess, 'localsingle'))

				h = figure(	'Units','points',...
							'Position',[0 0 800 800],...
							'Visible', 'off',...
							'Color', [1 1 1],...
							'Toolbar','none',...
							'NumberTitle','off',...
							'MenuBar','none' );

				movegui(h, 'center');
				set(h, 'Visible', 'on');

			else

				h = figure(	'Units','points',...
							'Position',[0 0 800 800],...
							'Visible', 'on',...
							'Color', [1 1 1],...
							'Toolbar','none',...
							'NumberTitle','off',...
							'MenuBar','none' );
			end	


			for subj = 1:nsub
				
				subject_data = [ ];
				
				for sess = aap.acq_details.selected_sessions

					temp = load(aas_getfiles_bystream(aap, subj, sess, 'metric_data'));
					metric_data = temp.metric_data;
					metric_data = metric_data.(metric_name);
					metric_data = metric_data(:); % force col vec
					subject_data = [ subject_data ; metric_data ];

				end
				
				metric_values = quantile(subject_data, quantiles_to_plot);
				plot(metric_values, 100*(1-quantiles_to_plot), 'b', 'LineWidth', 1);
				hold on;

				metric_subject_values(subj,:) = metric_values;
			
			end
			
			metric_mean_values = mean(metric_subject_values,1);
			
 			plot(metric_mean_values, 100*(1-quantiles_to_plot),'k', 'LineWidth', 4);
					
			% reference lines at 1%, 2% etc. TOI
						
			threshold_values = interp1(100*(1-quantiles_to_plot), metric_mean_values, TOI);

			% trim long tail
			a = axis; axis([a(1) 2*threshold_values(1) 0 100]);

			os = [ 1 1.5 1.5 1.5 1.5 1.5 1.5 1.5 ];	% text offsets
			
			for tindex = 1:length(TOI)
				
				a = axis;
				text_x = 0.9 * a(2);					
				text_y = TOI(tindex);

				hline = refline([0 text_y]);
				hline.Color = [0.8 0.8 0.8]; hline.LineWidth = 1;
				
				text(text_x, text_y+os(tindex), num2str(threshold_values(tindex),'%.3g'),'FontName', 'Helvetica', 'FontSize', 14, 'Color', [0.8 0.8 0.8]);
				
			end
					
			ylabel('data loss (%)');
			xlabel(metric_name);
			title([ metric_name ' vs. data loss (n = ' num2str(nsub) ')']);
			
			set(gca,'FontName','Helvetica','FontWeight','bold','FontSize',16,'LineWidth',2);
			
			if (nsub > 1)

				% boxplot the data losses across all subjects @ the mean threshold values

				boxdata = zeros(nsub,length(TOI));

				for subj = 1:nsub
					boxdata(subj,:) = interp1(metric_subject_values(subj,:), 100*(1-quantiles_to_plot), threshold_values);
				end

				axes('Position',[0.4 0.55 0.4 0.3]);
				bp = boxplot(fliplr(boxdata),'Labels',num2str(fliplr(threshold_values)','%0.3g'),'LabelOrientation','inline');
				a = axis; axis([a(1) a(2) 0 100 ]);
				set(bp,{'linew'},{2});
				set(findobj(gca,'Type','text'), 'FontSize', 14, 'HorizontalAlignment', 'center');
				ylabel('data loss (%)');
				title({['Data loss across all subjects at indicated threshold'];[' ']});
				set(gca,'FontName','Helvetica','FontWeight','bold','FontSize',12,'LineWidth',2);
				txt = findobj(gca,'Type','text');
				set(txt, 'FontWeight', 'bold');
				
			end

			fname = fullfile(aap.acq_details.root,[ 'threshold_' metric_name ]);
			set(h,'Renderer','opengl');
			set(findall(h,'Type','text'),'FontUnits','normalized');
			print(h, '-djpeg', '-r150', fname);
            close(h);

			% save threshold values for this metric for later stream output
						
			% fieldnames must match TOI = [ 1.0 2.0 5.0 10.0 20.0 40.0 ];
			
			metric_thresholds.(metric_name).onepercent		= threshold_values(1);
			metric_thresholds.(metric_name).twopercent		= threshold_values(2);
			metric_thresholds.(metric_name).fivepercent		= threshold_values(3);
			metric_thresholds.(metric_name).tenpercent		= threshold_values(4);
			metric_thresholds.(metric_name).twentypercent	= threshold_values(5);
			metric_thresholds.(metric_name).fortypercent	= threshold_values(6);
					
		end % loop over metric_name
		
		% desc metric metric_thresholds - this will overwrite
		% metric_data.mat with new struct containing threshold values

		outfile = fullfile(aas_getstudypath(aap), 'metric_thresholds.mat'); 
		save(outfile,'metric_thresholds');
		aap = aas_desc_outputs(aap, 'metric_thresholds', outfile);
 
		% done!

    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
		
end

end