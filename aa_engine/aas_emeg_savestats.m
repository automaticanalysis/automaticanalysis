function aas_emeg_savestats(D,newtits,newvals,newdefs,overwrite,fname,propagate)
% Helper function to save various statistics to (e.g. block-level) matfile
% aas_emeg_savestats(D,newtits,newvals,newdefs,overwrite,fname)
% Danny Mitchell 12/08
%
% Can have problems with parallel mode - file can get corrupted

global aaworker

if length(newtits)~=length(newvals) || length(newvals)~=length(newdefs)
    newtits, newvals, newdefs
    aas_log(aap,1,sprintf('\nFailed to save stats - list lengths do not agree.\n'))
end

if ~exist('fname','var'); fname='stats.mat'; end

if isstruct(D)
    [pth nam]=fileparts(D.fname);
    mout=fullfile(D.path,fname);
else
    [pth nam]=fileparts(D);
    mout=fullfile(pth,fname);
end
nam=regexprep(nam,'-','_');
lock=strrep(mout,'.mat','.lock');

if exist(mout,'file')
    numtries=0;
    while numtries<5000
        try 
            if exist('propagate','var') && ~isempty(aaworker)
                aas_propagatefrom(aaworker.master.hostname,mout);
                aas_propagatefrom(aaworker.master.hostname,lock);
            end
            if exist(lock,'file'); pause(0.1); gocatch % try for about 8 minutes?
            else
                diary(lock); diary off
                if exist('propagate','var') && ~isempty(aaworker); 
                    aas_propagateto(aaworker.master.hostname,lock);
                end
                try load(mout); 
                catch; delete(lock); gocatch % corrupt file?
                end
                break
            end
        catch
            numtries=numtries+1;
            if numtries==5000; debugnow; end
        end
    end
else stats.(nam).tits={};
end

if exist('stats','var') && isfield(stats,nam)
    %tits=tits'; vals=vals'; defs=defs'; % this fixed a teething problem - should now be redundant
    for t=1:length(newtits)
        switch sum(strcmp(stats.(nam).tits,newtits{t}))
            case 0; o=length(stats.(nam).tits)+1;
            case 1;
                if ~exist('overwrite','var') || ~overwrite
                    continue
                else
                    o=find(strcmp(stats.(nam).tits,newtits{t}));
                end
            otherwise
                aas_log(aap,1,sprintf('\nTwo values exist with same name.\n'))
        end
        stats.(nam).tits(o)=newtits(t);
        stats.(nam).vals(o)=newvals(t);
        stats.(nam).defs(o)=newdefs(t);
    end
else
    % convert row to column vectors if necessary
    stats.(nam).tits=newtits(:);
    stats.(nam).vals=newvals(:);
    stats.(nam).defs=newdefs(:);
end

numtries=0;
while numtries<100
    try 
        save(mout,'stats'); 
        if exist('propagate','var') && ~isempty(aaworker); 
            aas_propagateto(aaworker.master.hostname,mout); 
        end
        try % check it wrote ok, then remove lock
            if exist('propagate','var') && ~isempty(aaworker); 
                aas_propagatefrom(aaworker.master.hostname,mout); 
            end
            load(mout)
            if exist(lock,'file'); delete(lock); end
            break
        catch; gocatch % file may be corrupt? try again...
        end
    catch
        numtries=numtries+1;
        if numtries==100; debugnow; end
    end
end

return