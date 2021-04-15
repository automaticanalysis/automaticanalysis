% Main code for the multi-level reporting
% It distributes information between the different htmls.
% Tibor Auer MRC CBU Cambridge 2012-2013

function aap = aas_report_add(varargin)
aap = varargin{1};

if ~isfield(aap,'report')
    return
end

str = varargin{3};
if isempty(varargin{2})
    ptr = 'html_main';
else
    if isnumeric(varargin{2}) % Subject
        ptr = sprintf('html_S%02d',varargin{2});
        if ~isfield(aap.report,ptr)
            i = varargin{2};
            aap.report.(sprintf('html_S%02d',i)).fname = fullfile(aap.report.subdir,[aap.report.fbase sprintf('_S%02d.htm',i)]);
            aap = aas_report_add(aap,0,...
                sprintf('<a href="%s" target=_top>%s</a><br>',...
                aap.report.(sprintf('html_S%02d',i)).fname,...
                ['Subject: ' basename(aas_getsubjpath(aap,i))]));
            aap = aas_report_add(aap,i,['HEAD=Subject: ' basename(aas_getsubjpath(aap,i))]);            
        end
    else
        ptr = sprintf('html_%s',varargin{2});
    end
end
if isempty(str)
    aap.report.(ptr).fid = fopen(aap.report.(ptr).fname,'w');
elseif (numel(str) > 5) && strcmp(str(1:5), 'HEAD=')
    aap.report.(ptr).fid = fopen(aap.report.(ptr).fname,'w');
    
    % calculate path difference
    dpath = sum(strrep(spm_file(aap.report.(ptr).fname,'path'),spm_file(aap.report.html_main.fname,'path'),'')==filesep);
    dpath = reshape(char(arrayfun(@(x) '../',1:dpath,'UniformOutput',false))',1,[]);
    
    fprintf(aap.report.(ptr).fid,'<!DOCTYPE html>\n');
    fprintf(aap.report.(ptr).fid,'<html>\n');
    fprintf(aap.report.(ptr).fid,'<head><link rel="stylesheet" href="%saa_styles.css"></head>\n',dpath);
    if strcmp(ptr, 'html_main')
        fprintf(aap.report.(ptr).fid,'<meta charset="utf-8">\n');
    end
    fprintf(aap.report.(ptr).fid,'<body>\n');
    if strcmp(ptr, 'html_main')
        fprintf(aap.report.(ptr).fid,'<script src="https://d3js.org/d3.v5.min.js"></script>\n');
        fprintf(aap.report.(ptr).fid,'<script src="https://unpkg.com/@hpcc-js/wasm@0.3.11/dist/index.min.js"></script>\n');
        fprintf(aap.report.(ptr).fid,'<script src="https://unpkg.com/d3-graphviz@3.0.5/build/d3-graphviz.js"></script>\n');
    end
    fprintf(aap.report.(ptr).fid,'<table border=0>\n');
    fprintf(aap.report.(ptr).fid,'<td align=center width=100%%>\n');
    fprintf(aap.report.(ptr).fid,'%s\n',['<tr><td align=center><font size=+3><b>' str(6:end) '</b></font></tr>']);
    fprintf(aap.report.(ptr).fid,'</table>\n');
    fprintf(aap.report.(ptr).fid,'<a href="%s" target=_top>Main</a> &nbsp;-&nbsp;',aap.report.html_main.fname);
    fprintf(aap.report.(ptr).fid,'<a href="%s" target=_top>Subject list</a> &nbsp;-&nbsp;',aap.report.html_S00.fname);
    if isfield(aap.report,'html_moco')
        fprintf(aap.report.(ptr).fid,'<a href="%s" target=_top>Motion correction summary</a> &nbsp;-&nbsp;',aap.report.html_moco.fname);
    end
    if isfield(aap.report,'html_reg')
        fprintf(aap.report.(ptr).fid,'<a href="%s" target=_top>Registration summary</a> &nbsp;-&nbsp;',aap.report.html_reg.fname);
    end
    if isfield(aap.report,'html_C00')
        fprintf(aap.report.(ptr).fid,'<a href="%s" target=_top>First-level results</a>',aap.report.html_C00.fname);
    end
    if isfield(aap.report,'html_er')
        fprintf(aap.report.(ptr).fid,'<a href="%s" target=_top>MEEG epoch summary</a>',aap.report.html_er.fname);
    end
    fprintf(aap.report.(ptr).fid,'\n<hr class="rounded">\n');
    if strcmp(ptr, 'html_main')
        fprintf(aap.report.(ptr).fid,'<h2>Workflow</h2>\n');
        fprintf(aap.report.(ptr).fid,'<div id="workflow"></div>\n');
        fprintf(aap.report.(ptr).fid,'<script>\n');
        fprintf(aap.report.(ptr).fid,'d3.select("#workflow").graphviz()\n');
        fprintf(aap.report.(ptr).fid,'.renderDot(''digraph {''\n');
        
        % read dot data
        dotfid = fopen(spm_file(aap.report.(ptr).fname,'filename','aap_prov.dot'));
        lines = {};
        while ~feof(dotfid), lines{end+1} = fgetl(dotfid); end
        fclose(dotfid);
        lines([1 end]) = '';
        
        % write dot data to html
        for l = lines, fprintf(aap.report.(ptr).fid,['    +''' l{1}(2:end) '''\n']); end
        
        fprintf(aap.report.(ptr).fid,'    +''}'');\n');      
        fprintf(aap.report.(ptr).fid,'</script>\n');
    end
    fprintf(aap.report.(ptr).fid,'<table border=0>\n');
elseif strcmp(str, 'EOF')
    fprintf(aap.report.(ptr).fid,'</table>\n');
    fprintf(aap.report.(ptr).fid,'\n<hr class="rounded">\n');
    fprintf(aap.report.(ptr).fid,'</body>\n');
    fprintf(aap.report.(ptr).fid,'</html>\n');
    aap.report.(ptr).fid = fclose(aap.report.(ptr).fid);
else
    fprintf(aap.report.(ptr).fid,'%s\n',str);
    %     strcat(aap.report.(ptr).text,sprintf('%s\n',str));
end
end