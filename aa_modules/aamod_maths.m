function [aap,resp]=aamod_maths(aap,task,subj,sess)
resp='';

switch task
    case 'report'
        
    case 'doit'
        
        inputstreams=aas_getstreams(aap,'input'); istream = inputstreams{1};
        outputstreams=aas_getstreams(aap,'output'); ostream = outputstreams{1};
        
        % obtain data
        inputfnames = aas_getfiles_bystream(aap,aap.tasklist.currenttask.domain,[subj sess],istream);
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
        aap=aas_desc_outputs(aap,aap.tasklist.currenttask.domain,[subj sess],ostream,outputfnames);
        
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
