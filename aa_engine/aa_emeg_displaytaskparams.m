function pos=aa_emeg_displaytaskparams(taskstruc,varargin)

fields=fieldnames(taskstruc);

if isempty(varargin)
    [pth name]=fileparts(taskstruc.Path);
    spm_input(sprintf('Current options for task: %s',name),1,'d');
else
    spm_input(varargin{1},1,'d'); % use varargin as general heading
end

if length(fields)>1
    pos=2;
    for f=1:length(fields)
        value=taskstruc.(fields{f});
        if ischar(value);
        elseif isstruct(value)
            unfold2file(taskstruc.(fields{f}),fields{f});
            value=sprintf('Struct: %g member(s) & %g field(s).',size(value,2),length(fieldnames(value)));
        elseif iscellstr(value)
            temp=[char(value)'; blanks(size(value,2))];
            value=temp(:)';
        else
            try
                if iscell(value)
                    m=cell2mat(value);
                    value=reshape(m,length(value),[])';
                end
                value=mat2str(value);
            catch
            end
        end
        txt=sprintf('%s: %s',fields{f},value);
        % trim long values, because waiting for them to scroll is a bit annoying
        if length(txt)>53; txt=[txt(1:50) '...'];end
        spm_input(txt,pos,'d');
        pos=pos+1;
    end
end

return
