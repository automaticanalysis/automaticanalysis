function [aap]=aas_addbenchmark(aap)

callstack=dbstack;
callername=callstack(2).name;

try
    aap.internal.benchmarking.(callername).duration(end+1)=toc;
catch
    aap.internal.benchmarking.(callername).duration=toc;
end;