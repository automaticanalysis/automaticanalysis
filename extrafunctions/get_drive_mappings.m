function [mappings, error_msg] = get_drive_mappings()
% [mappings, error_msg] = get_drive_mappings()
% Return windows drive mappings as a nx2 cell array.
%
% Outputs:
% - mappings: nx2 cell array, with n>=0, and example output
%        {
%          'N:', '\\server1\dirX'
%          'Z:', '\\server2\dirY'
%        }
% - error_msg: empty string if successful

%% Developer note
% Getting the drive mappings is based on parsing the output of the net use
% command. This output depends on the Windows language settings. Especially,
% the Unavailable message may be a multiword text!
%
%     New connections will be remembered.
%
%
%     Status       Local     Remote                    Network
%
%     -------------------------------------------------------------------------------
%     OK           N:        \\server1\dir2            Microsoft Windows Network
%     Unavailable  S:        \\server1\archive         Microsoft Windows Network
%     OK           V:        \\server2\dir1            Microsoft Windows Network
%     OK           Z:        \\serverhere\mappingOfAQuiteLongPath\andSomeMore
%                                                     Microsoft Windows Network
%     The command completed successfully.

error_msg = '';
mappings = cell(0,2);

[status,result] = system('net use');

if status
    % Non-zero status means failure
    error_msg = result;
else
    % Split on newlines and find lines that have \\ in them
    lines = strsplit(result, '\n');
    for i = 1:length(lines)
        this_line = lines{i};
        if contains(this_line, '\\')
            % Split the line and find the part that has \\ in it
            line_parts = strsplit(this_line);
            idxRemote = find(contains(line_parts, '\\'), 1, 'first');
            idxDrive = idxRemote - 1;
            mappings(end+1,:) = line_parts(idxDrive:idxRemote); %#ok<AGROW> 
        end
    end
end

end