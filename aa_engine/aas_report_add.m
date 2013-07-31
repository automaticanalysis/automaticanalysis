% Main code for the multi-level reporting
% It distributes information between the different htmls.
% Tibor Auer MRC CBU Cambridge 2012-2013

function aap = aas_report_add(varargin)
aap = varargin{1};

if ~isfield(aap,'report')
    return
end

str = varargin{nargin};
if nargin == 2
    ptr = 'html_main';
else
    if isnumeric(varargin{2})
        ptr = sprintf('html_S%02d',varargin{2});
    else
        ptr = sprintf('html_%s',varargin{2});
    end
end
if isempty(str)
    aap.report.(ptr).fid = fopen(aap.report.(ptr).fname,'w');
elseif (numel(str) > 5) && strcmp(str(1:5), 'HEAD=')
    aap.report.(ptr).fid = fopen(aap.report.(ptr).fname,'w');
    fprintf(aap.report.(ptr).fid,'%s\n','<TABLE BORDER=0>');
    fprintf(aap.report.(ptr).fid,'%s\n','<td align=center width=100%>');
    fprintf(aap.report.(ptr).fid,'%s\n',['<tr><td align=center><font size=+3><b>' str(6:end) '</b></font></tr>']);
    fprintf(aap.report.(ptr).fid,'%s\n','</TABLE>');
    fprintf(aap.report.(ptr).fid,'%s\n','<TABLE BORDER=0>');
    fprintf(aap.report.(ptr).fid,'%s',...
        sprintf('<a href="%s" target=_top>Main</a> &nbsp;-&nbsp',aap.report.html_main.fname));
    fprintf(aap.report.(ptr).fid,'%s',...
        sprintf('<a href="%s" target=_top>Subject list</a> &nbsp;-&nbsp;',aap.report.html_S00.fname));
    fprintf(aap.report.(ptr).fid,'%s',...
        sprintf('<a href="%s" target=_top>Motion correction summary</a> &nbsp;-&nbsp;',aap.report.html_moco.fname));
    fprintf(aap.report.(ptr).fid,'%s',...
        sprintf('<a href="%s" target=_top>Registration summary</a> &nbsp;-&nbsp;',aap.report.html_reg.fname));
    fprintf(aap.report.(ptr).fid,'%s',...
        sprintf('<a href="%s" target=_top>First-level contrasts</a>;',aap.report.html_C00.fname));
    fprintf(aap.report.(ptr).fid,'%s\n','</TABLE><hr>');
    fprintf(aap.report.(ptr).fid,'%s\n','<TABLE BORDER=0>');    
elseif strcmp(str, 'EOF')
    fprintf(aap.report.(ptr).fid,'%s\n','</TABLE><hr>');    
    aap.report.(ptr).fid = fclose(aap.report.(ptr).fid);
else
    fprintf(aap.report.(ptr).fid,'%s\n',str);
    %     strcat(aap.report.(ptr).text,sprintf('%s\n',str));
end
end