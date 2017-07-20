% Adds an event to a model
% function aap=aid_addcontrast(aap, modulename, subject, contrastname, matrix)
% 
% modulename = name of module (e.g.,'aamod_MVPaa_roi_1st') for contrast
% subject = subject for this contrast
% contrastname = name of the contrast matrix
% matrix = contrast matrix itself

function aap=mvpaa_addcontrast(aap,modulename,subject,contrastname,matrix)

% Get number from end of module name if present in format _%05d (e.g, _00001)
if (length(modulename>6)) %#ok<ISMT>
    moduleindex=str2num(modulename(end-4:end));
    if (~strcmp(['_' sprintf('%05d',moduleindex)],modulename(length(modulename)-5:end)))
        moduleindex=1;
    else
        modulename=modulename(1:length(modulename)-6);
    end
else
    moduleindex=1;
end

% find model that corresponds and add event to this if it exists
whichmodel=strcmp({aap.tasksettings.(modulename)(moduleindex).model.subject},subject);

if (~any(whichmodel))
    emptymod=[];
    emptymod.subject=subject;
    emptymod.contrast.name=contrastname;
    emptymod.contrast.matrix=matrix;    
    aap.tasksettings.(modulename)(moduleindex).model(end+1)=emptymod;
else
    aap.tasksettings.(modulename)(moduleindex).model(whichmodel).contrast(end+1).name=contrastname;
    aap.tasksettings.(modulename)(moduleindex).model(whichmodel).contrast(end).matrix=matrix;
end