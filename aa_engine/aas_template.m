function aap=aas_template(aap,templatefn,parameters,outputfn)

fid=fopen(templatefn);
fid_out=fopen(outputfn,'w');

while (~feof(fid))
    lne=fgetl(fid);
    
    pos_start=strfind(lne,'{{');
    pos_end=strfind(lne,'}}');
    if (length(pos_start)~=length(pos_end))
        aas_log(aap,true,sprintf('Unmatching {{ and }} in line %s of template file %s',lne,templatefn));
    end;
    outlne='';
    firstind=1;
    for curlyind=1:length(pos_start)
        parmname=strtrim(lne((pos_start(curlyind)+2):(pos_end(curlyind)-2)));
        fprintf('parmname %s value %f\n',parmname, parameters.(parmname));
        if (pos_start(curlyind)>1)
            outlne=[outlne lne(firstind:(pos_start(curlyind)-1))];
        end;
        value=parameters.(parmname);   
        if isnumeric(value)
            value=num2str(value);
        end;
        outlne=[outlne value];
        firstind=pos_end(curlyind)+2;
    end;
    if (firstind<=length(lne))
        outlne=[outlne lne(firstind:end)];
    end;
       
    fprintf(fid_out,'%s\n',outlne);
    
end;

fclose(fid);
fclose(fid_out);