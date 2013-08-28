% Convert the whole report between linux and windows (MRC-CBSU!)
% Tibor Auer MRC CBU Cambridge 2012-2013

function aas_report_convert(studyroot,l2w)
% l2w: Linux2Win
%     - true: Linux to Win - defualt
%     - false: Win to Linux

patt = {'/imaging',... % Linux path
    '/////cbsu/data/imaging'... % Win path
    }; 

try, cd(studyroot); catch; end
try, l2w; catch; l2w = true; end

p1 = patt{(~l2w)+1};
p2 = patt{l2w+1};

% First, load AAP structure
load('aap_parameters_reported');

fields = fieldnames(aap.report);
for f = 1:numel(fields)
    fname = fields{f};
    if isempty(strfind(fname,'html_'))
        continue;
    end
    fname = aap.report.(fname).fname;
    fid = fopen(fname,'r');
    buff = {};
    while ~feof(fid)
        buff{end+1} = fgetl(fid);
    end
    fclose(fid);
    fid = fopen(fname,'w');
    for l = 1:numel(buff)
        if (l2w && isempty(strfind(buff{l},patt{2}))) || (~l2w && ~isempty(strfind(buff{l},patt{2}))),
            buff{l} = strrep(buff{l},p1,p2);
        end
        fprintf(fid,'%s\n',buff{l});
    end
    fclose(fid);    
end