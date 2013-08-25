function [aap outstruct]= aas_xmltostruct(aap,inxml)

xmlfn=aas_gettempfilename();
fid=fopen(xmlfn,'w');
fprintf(fid,'%s',strtrim(inxml));
fclose(fid);
outstruct=xml_read(xmlfn);
delete(xmlfn);

