function aas_time_elapsed
% TIME_ELAPSED shows time past from tic...

% minimal aap structure foe aas_log
aap.options.verbose = 2;
aap.options.email = '';
aap.gui_controls.usecolouroutput = true;

time = toc; hours = floor(time/3600); mins = floor(rem(time,3600)/60); secs = floor(rem(rem(time,3600),60));
aas_log(aap,false,sprintf('\n\t \t \t Elapsed time = %dh:%dm:%ds',hours, mins, secs),'blue')
c = clock;
aas_log(aap,false,sprintf('\t \t \t ...at = %02.0f:%02.0f:%02.0f (%d/%d/%d)',c(4), c(5), c(6), c(3), c(2), c(1)),'blue')

