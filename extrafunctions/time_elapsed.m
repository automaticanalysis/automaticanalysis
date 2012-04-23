function time_elapsed
% TIME_ELAPSED shows time past from tic...

time = toc; hours = floor(time/3600); mins = floor(rem(time,3600)/60); secs = floor(rem(rem(time,3600),60));
try
    cprintf('blue', '\n\t \t \t Elapsed time = %dh:%dm:%ds \n',hours, mins, secs)
catch
fprintf('\n\t \t \t Elapsed time = %dh:%dm:%ds \n',hours, mins, secs)
end
c = clock;
try
    cprintf('blue', '\t \t \t ...at = %02.0f:%02.0f:%02.0f (%d/%d/%d) \n',c(4), c(5), c(6), c(3), c(2), c(1))
catch
    fprintf('\t \t \t ...at = %02.0f:%02.0f:%02.0f (%d/%d/%d) \n',c(4), c(5), c(6), c(3), c(2), c(1))
end

