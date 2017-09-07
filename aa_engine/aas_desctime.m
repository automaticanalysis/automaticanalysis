function desc=aas_desctime(secs)
desc=[];
% possibly a little pessimistic, and doesn't fit on a sensibly sized
% progress bar
% millenia=floor(secs/24/3600/365/1000);
% if (millenia>0)
%     desc=[desc sprintf('%d millenia ',millenia)];
%     secs=secs-millenia*24*3600*365*1000;
% end;
% decades=floor(secs/24/3600/365/10);
% if (decades>0)
%     desc=[desc sprintf('%d decades ',decades)];
%     secs=secs-decades*24*3600*365*10;
% end;
years=floor(secs/24/3600/365);
if (years>0)
    desc=[desc sprintf('%d years ',years)];
    secs=secs-years*24*3600*365;
end;
weeks=floor(secs/24/3600/7);
if (weeks>0)
    desc=[desc sprintf('%d weeks ',weeks)];
    secs=secs-weeks*24*3600*7;
end;
days=floor(secs/24/3600);
if (days>0)
    desc=[desc sprintf('%d days ',days)];
    secs=secs-days*24*3600;
end;
hours=floor(secs/3600);
secs=secs-hours*3600;
mins=floor(secs/60);
secs=round(secs-mins*60);
desc=[desc sprintf('%d:%02d:%02d',hours,mins,secs)];
end