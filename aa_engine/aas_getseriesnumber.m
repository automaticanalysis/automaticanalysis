function [aap seriesnum]=aas_getseriesnumber(aap,fn)

switch (aap.directory_conventions.seriesnamingconvention)
    case 'CABI'
        findunderscore=find(fn=='_');
        seriesnum=str2num(fn((findunderscore(end)+1):end));
    case 'CBU'
        seriesnum=str2num(fn(8:10));
    case 'AWS'
        seriesnum=str2num(fn(8:11));
end;