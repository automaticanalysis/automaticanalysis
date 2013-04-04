function cg_vbm8_debug
%cg_vbm8_debug	print debug information for SPM8 and VBM8
%
% FORMAT cg_vbm8_debug
%
%__________________________________________________________________________
% Christian Gaser
% $Id: cg_vbm8_debug.m 404 2011-04-11 10:03:40Z gaser $

rev = '$Rev: 404 $';

% print last error
fprintf('\nLast error message:\n');
fprintf('-------------------------------------------------------------------------------------\n');
fprintf('-------------------------------------------------------------------------------------\n');
try
	er = lasterror;
	fprintf('%s\n',er.message);
	if isfield(er,'stack')
		for i=1:length(er.stack)
			fprintf('%s at line %g\n',char(er.stack(i).file),er.stack(i).line);
		end
	end
catch
	fprintf('%s\n',lasterr);
end

fprintf('-------------------------------------------------------------------------------------\n');
fprintf('-------------------------------------------------------------------------------------\n');

fprintf('\nVersion information:\n');
fprintf('-------------------------------------------------------------------------------------\n');

ver

