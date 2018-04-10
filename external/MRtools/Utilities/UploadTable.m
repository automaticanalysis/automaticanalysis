function UploadTable(dat,destTable)
% UploadTable(dat,'aschultz.BlindPenData');

cnames = dat.Properties.VariableNames;

R = ['Create Table ' destTable ' ('];
I = ['INSERT INTO ' destTable ' ('];
I2 = 'VALUES(';
for ii = 1:numel(cnames)
    I = [I cnames{ii} ','];
    
    if ~isempty(findstr('Date',cnames{ii}))
       R =  [R cnames{ii} ' Date,'];
       I2 = [I2 '"{S}",'];
       continue
    end
    
    if ~isempty(findstr('Session_ID',cnames{ii}))
       R =  [R cnames{ii} ' INT(11) NOT NULL Primary Key,'];
       I2 = [I2 '"{S}",'];
       continue
    end
    
    if isnumeric(dat.(cnames{ii}))
       R =  [R cnames{ii} ' float,'];
       I2 = [I2 '"{Sn}",'];
       continue
    end
    
    if iscell(dat.(cnames{ii}))
       R =  [R cnames{ii} ' varchar(255),'];
       I2 = [I2 '"{S}",'];
       continue
    end
    
    if isobject(dat.(cnames{ii}))
       R =  [R cnames{ii} ' varchar(255),'];
       dat.(cnames{ii}) = cellstr(dat.(cnames{ii}));
       I2 = [I2 '"{S}",'];
       continue
    end
end

% I = [I(1:end-1) ') ' I2(1:end-1) ')'];
I = [I(1:end-1) ') VALUES'];
R = [R(1:end-1) ')'];

II = I;
for ii = 1:size(dat,1)
    II = [II '('];
    for jj = 1:size(dat,2)
        tmp = dat.(cnames{jj})(ii);
        
        if iscell(tmp)
            II = [II '"' tmp{1} '",'];
            continue
        end
        
        if isnumeric(tmp)
            II = [II '"' num2str(tmp) '",'];
            continue
        end
        
        if isempty(tmp)
            II = [II 'NULL,'];
            continue
        end
        
        disp('got here.  Need more config'); keyboard;
    end
    
    II = [II(1:end-1) '), '];
end

II = II(1:end-2);

% INSERT INTO tbl_name (a,b,c) VALUES(1,2,3),(4,5,6),(7,8,9);
%%%
DataCentral(['DROP TABLE IF EXISTS ' destTable ]);
DataCentral(R);
DataCentral(II);
%%%



