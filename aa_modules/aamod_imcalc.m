% AA module - use imcalc to create a calculated image
%
% The module will run each function included in 
% aap.options.imcalc.functions separately for each subject
% and session.
%
% All files starting with the prefix in aap.options.imcalc.filefilt
% will be included, and the output will be an image with 'aa_imcalc'
% and the function name appended to the start of the first included
% image name.
%
% e.g.
% aap.options.imcalc.functions = {'mean(X)','var(X)'}
% aap.options.imcalc.filefilt = 'f'
% will create a mean image aa_imcalc_mean_f* and a variance
% image aa_imcalc_var_f*
%
% aap.options.imcalc can also include a field 'flags' with
% 4 values indicating:
%   whether or not to read image values to data matrix (0/1)
%   whether to implicitly mask zeros (0/1)
%   data type of output image (see spm_slice_vol)
%   interpolation hold
% The first value must be set as 1 for this module to work. If
% the field is not present aap.options.imcalc.flags defaults to
% {1;0;4;1} - read values into matrix, don't mask zeros, 16 bit
% integer, trilinear interpolation
%
% i=subject num
% j=session num

function [aap,resp]=aamod_meanvarimg(aap,task,i,j)

resp='';

switch task
    case 'domain'
        resp='session';   % this module needs to be run once per session
    case 'whentorun'
        resp='justonce';  % should this be run everytime or justonce?
    case 'description'
        resp = 'Imcalc ';
        for c = 1:size(aap.options.imcalc.functions,2)
            fn = aap.options.imcalc.functions{c};
            fn =fn(1:findstr('(X)',fn)-1);
            if c>1; resp = [resp ', '];end
            resp = [resp fn];
        end
    case 'summary'
        resp='Create '
        for c = 1:size(aap.options.imcalc.functions,2)
            fn = aap.options.imcalc.functions{c};
            fn =fn(1:findstr('(X)',fn)-1);
            if c>1; resp = [resp ', '];end
            resp = [resp fn];
        end

        resp = [resp 'images using imcalc\n'];
    case 'report'
        aap.report.html=strcat(aap.report.html,'<table><tr><td>');
        aap=aas_report_addimage(aap,fullfile(aas_getsesspath(aap,i,j),'diagnostic_aamod_imcalc(function).jpg'));
        aap.report.html=strcat(aap.report.html,'</td></tr></table>');
    case 'doit'
        subjpath=aas_getsubjpath(aap,i);
        sesspath=aas_getsesspath(aap,i,j);

        if (~length(dir(sesspath)))
            aas_log(aap,1,sprintf('Problem finding directory for session\n%s',sesspath));
        end;

        dirn = aas_getsesspath(aap,i,j);
        imgs=aas_getimages(aap,i,j,aap.options.imcalc.filefilt,0,inf);
        [p f e] =fileparts(imgs(1,:));
        
        if isfield(aap.options.imcalc,'flags')
            if isempty(aap.options.imcalc.flags)
                aap.options.imcalc.flags = {1;0;4;1};
            else
                aap.options.imcalc.flags{1} = 1;
            end
        else
            aap.options.imcalc.flags = {1;0;4;1};
        end

        for c = 1:size(aap.options.imcalc.functions,2)
            fn = aap.options.imcalc.functions{c};
            outfile =fullfile(dirn, ['aa_imcalc_' fn(1:findstr('(X)',fn)-1) '_' f e]);

            spm_imcalc_ui(imgs,outfile,fn,aap.options.imcalc.flags);
        end

    case 'checkrequirements'

    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end;