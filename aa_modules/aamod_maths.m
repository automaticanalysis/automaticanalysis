function [aap,resp]=aamod_maths(aap,task,varargin)
resp='';

switch task
    case 'report'
        
    case 'doit'
        
        inputstreams=aas_getstreams(aap,'input'); istream = inputstreams{1};
        outputstreams=aas_getstreams(aap,'output'); ostream = outputstreams{1};
        
        % obtain data
        inputfnameslist = cellstr(aas_getfiles_bystream(aap,aap.tasklist.currenttask.domain,cell2mat(varargin),istream));
        if ~isempty(aas_getsetting(aap,'inputselection'))
            inputfnameslist = inputfnameslist(cellfun(@(i) ~isempty(regexp(i,aas_getsetting(aap,'inputselection'),'match')), inputfnameslist));
        end
        
        % run
        outputfnameslist = '';
        for inp = 1:numel(inputfnameslist)
            inputfnames = inputfnameslist{inp};
            V = spm_vol(inputfnames);
            dat = spm_read_vols(V);
            
            % process operation string
            ops = aap.tasklist.currenttask.settings.operation;
            if ~iscell(ops), ops = {ops}; end
            
            for o = 1:numel(ops)
                op = ops{o};
                
                op(op==' ') = '';
                par_open = find(op=='(');
                par_close = find(op==')');
                
                funcstr = op(1:par_open-1);
                
                func = str2func(funcstr);
                
                args = textscan(op(par_open+1:par_close-1),'%s','Delimiter',','); args = args{1};
                
                % substitute data
                args{strcmp(args,'X')} = dat;
                
                % special cases
                args(strcmp(args,'[]')) = {[]};
                
                for a = 1:numel(args)
                    if ischar(args{a}) && ~isempty(str2num(args{a})), args{a} = str2num(args{a}); end
                end
                
                dat = func(args{:});
            end
            % Describe outputs
            outputfnames = spm_file(inputfnames,'prefix',[funcstr '_']);
            V = V(1);
            V.dt = spm_type('float32');
            nifti_write(outputfnames,dat,funcstr,V);
            outputfnameslist(inp,:) = outputfnames;
        end
        aap=aas_desc_outputs(aap,aap.tasklist.currenttask.domain,cell2mat(varargin),ostream,outputfnameslist);
        
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end
end

%%%%%%%%%%%%%%%%%%%%% MATHS %%%%%%%%%%%%%%%%%%

function o = div(x,c)
o = x/c;
end

function o = mul(x,c)
o = x*c;
end

function o = add(x,c)
o = x+c;
end

function o = sub(x,c)
o = x-c;
end

% binarise
function o = bin(x)
x(isnan(x)) = 0;
o = logical(x);
end

% split a multi-label volume (i.e. one volume with multiple non-overlapping labels) into singe-label volumes
function o = splitbylabels(x)
o = {};
for l = unique(x(x>0))'
    o{end+1} = x == l;
end
o = cat(ndims(o{1})+1,o{:});
end

% identical volume (usefull for input selection)
function o = ident(x)
o = x;
end