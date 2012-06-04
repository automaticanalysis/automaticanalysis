function [aap resp]=aa_feval_withindices(mfile_alias,aap,task,indices)
switch(length(indices))
    case 0
        [aap,resp]=aa_feval(mfile_alias,aap,task);
    case 1
        [aap,resp]=aa_feval(mfile_alias,aap,task,indices(1));
    case 2
        [aap,resp]=aa_feval(mfile_alias,aap,task,indices(1),indices(2));
    case 3
        [aap,resp]=aa_feval(mfile_alias,aap,task,indices(1),indices(2),indices(3));
    case 4
        [aap,resp]=aa_feval(mfile_alias,aap,task,indices(1),indices(2),indices(3),indices(4));
    case 5
        [aap,resp]=aa_feval(mfile_alias,aap,task,indices(1),indices(2),indices(3),indices(4),indices(5));
    case 6
        [aap,resp]=aa_feval(mfile_alias,aap,task,indices(1),indices(2),indices(3),indices(4),indices(5),indices(6));
    otherwise
        aas_log(aap,true,'More than 6 indices not currently supported.');
end;
end