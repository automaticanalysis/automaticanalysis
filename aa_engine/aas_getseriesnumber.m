function [aap seriesnum]=aas_getseriesnumber(aap,fn)

format = strrep(regexprep(aap.directory_conventions.seriesoutputformat,'%.*d','[0-9]*'),'**','*');

seriesnum = regexp(fn,format,'match');
seriesnum = seriesnum{end};

seriesnum = regexp(seriesnum,'[0-9]*','match');
seriesnum = str2double(seriesnum{1});
