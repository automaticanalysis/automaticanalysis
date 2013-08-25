% Automatic analysis - pvshow.pl
% This function is no longer used.
% It parse's output of Bruker dump tool pvshow.pl by Matthew Brett
% Rhodri Cusack MRC CBU Cambridge 2004

function  [aap,seriesnum,acqtype,acqtime,reps,recoed]=aas_parsepvshow(aap,rawdatasetpath)

cmd=sprintf('pvshow.pl %s',rawdatasetpath);
[s,w]=aas_shell(cmd);

if (s | strmatch('Found no matching data on searchpath',w))
    aas_log(aap,1,sprintf('Error running pvshow.pl\ncommand: %s\ngave error: %s',cmd,w));   
end;

lineends=findstr(10,w);

delims=['::;;;;']; 
lines={};
linestart=1;

seriesnum=[];
acqtype=[];
acqtime=[];
reps=[];
recoed=[];

for i=1:length(lineends)
    thisline=w(linestart:(lineends(i)-1));
    lines{i}={};
    rem=[thisline ';'];
    % check this line contains the path - otherwise it is probably an
    %   error and should be ignored
    [pth fn ext]=fileparts(rawdatasetpath);
    if (length(findstr(thisline,pth))<1)
        aas_log(aap,0,sprintf('*** WARNING strange line in pvshow.pl output "%s" - expected to find %s in it\n',thisline,pth))
    else
        
        for j=1:length(delims)
            ind=findstr(delims(j),rem);
            if (length(ind)<1) 
                ass_log(aap,1,sprintf('Error parsing pvshow.pl output. Has a new version been installed? If so, I dont understand it!'));
            end;
            myfield=rem(1:(ind(1)-1));
            rem=rem((ind(1)+1):length(rem));
            switch j
                case 2
                    seriesnum(i)=str2num(myfield);
                case 3
                    acqtype{i}=myfield;
                case 4
                    acqtime{i}=myfield;
                case 5
                    reps(i)=sscanf(myfield(7:length(myfield)),'%d');
                case 6
                    if (findstr(myfield,'reco: Yes'))
                        recoed(i)=1;
                    else
                        recoed(i)=0;
                    end;
            end;
        end; 
        
    end;        
    linestart=lineends(i);    
end;
