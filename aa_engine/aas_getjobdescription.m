function [resp]=aas_getjobdescription(aap,title,joblist,joberrors)

maxlist=40;

resp=['<table><tr><td>' num2str(length(joblist)) ' ' title '</td></tr>'];
if (~isempty(joblist))
    for i=1:min(maxlist,length(joblist))
        myerror='';
        for j=1:length(joberrors{joblist(i)})
            myerror=[myerror sprintf('<a href="joberrors/J%d_%d.html">e</a>&nbsp;', joblist(i),j)];
        end;
        shortname=aap.internal.taskqueue(joblist(i)).stagename;
        shortname=strrep(shortname,'aamod_','');
        shortname=strrep(shortname,'emeg_','');        
        resp=[resp sprintf('<tr><td><a href="#" onclick="alert(''%s'');" title="%s">J%d-%s</a> <i>%s</i></td></tr>',aap.internal.taskqueue(joblist(i)).description,aap.internal.taskqueue(joblist(i)).description,joblist(i),shortname,myerror)];
    end;
            
    if (length(joblist)>maxlist)
        resp=[resp sprintf('<tr><td>[%d more]...</td></tr>',length(joblist)-maxlist)];
    end;
end;
resp=[resp '</table>'];        