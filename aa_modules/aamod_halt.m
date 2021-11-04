function [aap, resp] = aamod_halt(aap,task)
%
% halt execution
%
% CHANGE HISTORY
%
% 01/2019 [MSJ] -- new
%

resp = '';

switch task
	
    case 'checkrequirements'
    case 'report'
    case 'doit'
		
        aas_log(aap, true, sprintf('INFO (%s): Halting execution as requested. This is not an error.', mfilename));
		        
    otherwise
        
        aas_log(aap, true, sprintf('%s:Unknown task %s', mfilename, task));
		
end

end
