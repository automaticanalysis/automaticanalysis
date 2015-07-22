function gen_object_display( obj_struct,indent )
%
% gen_object_display - general function to display an object's content
%
% format:   gen_object_display( obj_struct,indent )
%
% input:    obj_struct  - a copy of the object stored inside a structure
%           indent      - amount of "indent" when printing to the screen
%
% output:   to the screen
%
% example:  gen_object_display( struct( my_object_handle) );
%           gen_object_display( ny_structure );
%
% Correction History:
%  2006-11-01 - Jarek Tuszynski - added support for struct arrays

%% handle insufficient input
if ( nargin == 0 )
    help gen_object_display;
    return;
elseif (nargin == 1)
    indent = 1;
end

%% check input for errors
% if ~isstruct( obj_struct )
%     fprintf( '\n\n\tMake sure that ''obj_struct'' is a struct type\n' );
%     return
% end

% if (iscell( obj_struct ))
%   for i =1:length(obj_struct)
%     gen_object_display( obj_struct{i},indent + 2 );
%   end
%   return
% end
if ~isstruct( obj_struct )
  space = sprintf( sprintf( '%%%ds',indent ),' ' );
  fprintf( '   %s', space);
  disp(obj_struct);
  return
end

% find the longest name
field_list      = fieldnames( obj_struct );
max_strlen      = 0;
for idx = 1:length( field_list )
    max_strlen  = max( max_strlen,length(field_list{idx}) );
end

%% setup the display format (spacing)
space       = sprintf( sprintf( '%%%ds',indent ),' ' );
name_format = sprintf( '   %s%%%ds: ', space, max_strlen );
max_displen = 110 - max_strlen - indent;

%% display each field, if it is not too long
for iItem = 1:length( obj_struct ) % loop added  by JT
  for idx = 1:length( field_list )

    % prepare field name to be displayed
    name = sprintf( name_format,field_list{idx} );
    %temp = getfield( obj_struct,field_list{idx} ); % original by OG
    temp = obj_struct(iItem).(field_list{idx});    % modification by JT

    % proceed according the variable's type
    switch (1)
    case islogical( temp ), % case added by JT
        if (temp)
            fprintf( '%strue\n',name );
        else
            fprintf( '%sfalse\n',name );
        end
    case ischar( temp ),
        if (length(temp)<max_displen )
            fprintf( '%s''%s''\n',name,temp' );
        else
            fprintf( '%s[%dx%d char]\n',name,size(temp,1),size(temp,2) );
        end
    case isnumeric( temp ),
        if (size( temp,1 )==1 )
            temp_b = num2str( temp );
            if (length(temp_b)<max_displen )
                fprintf( '%s[%s]\n',name,temp_b );
            else
                fprintf( '%s[%dx%d double]\n',name,size(temp,1),size(temp,2) );
            end
        else
            fprintf( '%s[%dx%d double]\n',name,size(temp,1),size(temp,2) );            
        end
    case iscell( temp ),    
      if (length(temp)<10)
        fprintf( '%s[%dx%d cell] = \n',name,size(temp,1),size(temp,2) );
        %disp(temp)
        for i =1:length(temp)
          gen_object_display( temp{i},indent + max_strlen + 2 );
          fprintf('\n');
        end
      else
        fprintf( '%s[%dx%d cell]\n',name,size(temp,1),size(temp,2) );
      end
    case isstruct( temp ),  fprintf( '%s[%dx%d struct]\n',name,size(temp,1),size(temp,2) );
        if (indent<80)
            gen_object_display( temp,indent + max_strlen + 2 );
        end
    case isobject( temp ),  fprintf( '%s[inherent object]\n',name );
        if (indent<80)
            cmd = sprintf( 'display( obj_struct.%s,%d );',field_list{idx},indent + max_strlen + 2 );
            eval( cmd );
        end
    otherwise,
        fprintf( '%s',name );
        try
            fprintf( temp );
        catch
            fprintf( '[No method to display type]' );
        end
        fprintf( '\n' );
    end
  end
  if (length(obj_struct)>1), fprintf('\n'); end % added by JT
end                                             % added by JT