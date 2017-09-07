function [aap, resp]=aa_feval_withindices(mfile_alias,aap,task,indices)

ci = num2cell(indices);
[aap,resp]=aa_feval(mfile_alias,aap,task,ci{:});
end