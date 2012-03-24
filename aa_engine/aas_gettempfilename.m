function tmpname=aas_gettempfilename()
while(1)
    tmpname=tempname;
    if (~exist(tmpname,'file'))
        break;
    end;
end;