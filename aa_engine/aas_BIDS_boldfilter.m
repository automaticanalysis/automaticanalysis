function aap = aas_BIDS_boldfilter(aap,boldfilter)

% remove subject sessions from the the current analysis that
% fail a specified filter based on json header information
% in the BIDS /func directory
%
% this function modifies aap.acq_details.subjects.seriesnumbers
%
% NB: call this *AFTER* aas_processBIDS
%
% --- boldfilter syntax ---
% 
%       struct with fields:
% 
%         fieldname: 'meanSNR'
%                op: '>'        ('>','<', or '=')
%             value: 1.23
% 
%      to define one filter appled to all subjects, or
% 
%       struct with fields:
% 
%         fieldname: 'meanSNR'
%                op: '>'
%             sub01: 1.23
%             sub02: 4.56 (etc...)
% 
%      to define a subject-specific filter.
% 
%      here we have assumed "meanSNR" is a field that appears in the json
% 
%      "op" can be '<','>', or '=' ('==' is also recognized)
% 
%      additionally, use the op '><' to test the value against a range:
% 
%       struct with fields:
% 
%         fieldname: 'meanSNR'
%                op: '><'
%             value: [1.2 3.4]
% 
%       or struct with fields:
% 
%         fieldname: 'meanSNR'
%                op: '><'
%             sub01: [1.2 3.4]
%             sub02: [5.6 7.8] (etc...)
% 
% 
%     NB: for subject-specific filtering, the special characters in the
%     BIDS subject identifier do not appear (e.g. use "sub01" not "sub-01")
%     This is because struct fieldnames can only contain alphanumerics
% 
%     NB: equality ("=") should not be used with floating point values
%     because of roundoff uncertainty. Use >< instead.
%

% NEW!

%         fieldname: 'meanSNR'
%                op: 'M'        ('use session with Maximum meanSNR')

%         fieldname: 'meanSNR'
%                op: 'm'        ('use session with minimum meanSNR')



subject_names = { aap.acq_details.subjects.subjname };


if (strcmp(boldfilter.op,'M') || strcmp(boldfilter.op,'m'))
    
    % best and worst require a screening pass where the best (or worst)
    % session is identified and the filter is modified to pick it out
    % in the main pass
    
    for subject_index = 1:numel(subject_names)

        subject_values = [];

        for session_index = 1:numel(aap.acq_details.subjects(subject_index).seriesnumbers)

            json_fname = aap.acq_details.subjects(subject_index).seriesnumbers{session_index}{1}.hdr;

            header = loadjson(json_fname);

            if ~isfield(header, boldfilter.fieldname)

                aas_log(aap,true,sprintf('WARNING: specified BIDSboldfilter field %s is not defined in *_bold.json. Halting.', boldfilter.fieldname));

            else

              subject_values(end+1) = header.(boldfilter.fieldname);

            end

        end

        temp = subject_names{subject_index};
        temp = split(temp,'_');
        stripped_subjname = strrep(temp{1},'-',''); 

        if strcmp(boldfilter.op,'M'); target_value = max(subject_values); end
        if strcmp(boldfilter.op,'m'); target_value = min(subject_values); end
        
        % define a delta to bracket the selection
        
        delta = 1;
        if (length(subject_values) > 1)
            delta = min(abs(diff(subject_values)))/2;
        end

        boldfilter.(stripped_subjname) = [ target_value-delta target_value+delta ];

    end
    
    % we can now get the session of interest using a range op
    
    boldfilter.op = '<>';

end


      

subject_kill_list = []; % subjects with no sessions that pass filter

% aap.acq_details.input.combinemultiple is a bit tricky
%
% if combinemultiple == false, subject_names will be sub-01_ses-01
% sub-01_ses-02 ... sub_02_ses01 sub02_ses01 ...
%
% if combinemultiple == true, subject names will be sub-01, sub-02 ...
% and each will have multiple entries in seriesnumber

for subject_index = 1:numel(subject_names)
    
    keepers = []; % acq_details.subject(subj).seriesnumbers entries that pass the filter
        
    % if combinemultiple == false, each subject should have one session
    
    for session_index = 1:numel(aap.acq_details.subjects(subject_index).seriesnumbers)

        json_fname = aap.acq_details.subjects(subject_index).seriesnumbers{session_index}{1}.hdr;

        header = loadjson(json_fname);

        if ~isfield(header, boldfilter.fieldname)

            aas_log(aap,true,sprintf('WARNING: specified BIDSboldfilter field %s is not defined in *_bold.json. Halting.', boldfilter.fieldname));

        else

            % apply filter to this file
            
            % '-' can't appear in struct fieldnames used in the boldfilter
            % also, for combinemultiple = false, the subject_name will have
            % the session ID appended, which isn't used in the boldfilter
            
            temp = subject_names{subject_index};
            temp = split(temp,'_');
            stripped_subjname = strrep(temp{1},'-',''); 

            value = header.(boldfilter.fieldname);

            op = boldfilter.op;

            if strcmp(op,'==');op='=';end

            if (isfield(boldfilter,stripped_subjname))
                test_value = boldfilter.(stripped_subjname);
            else
                test_value = boldfilter.value;
            end

            if (length(test_value) == 1)
                pass = eval_BIDSboldfilter(aap,value,op,test_value,[]);
            else
                op = '><';
                pass = eval_BIDSboldfilter(aap,value,op,test_value(1),test_value(2));
            end

            % if this epi passes the filter, gets to stay
            % in acq_details.subjects(subj) (see below)
            
            % if combinemultiple == true, pass/fail applies to series numbers
            % if combinemultiple == false, pass/fail applies to subject
            % (bc each subject/session is treated like a separate subject)

            if (pass==true); keepers(end+1) = session_index; end

        end % isfield(boldfilter)

    end % loop over series numbers
        
    % if no epi passed the filter, then keepers will be empty and we
    % should delete the subject 
    
    % if keepers is not empty, delete sessions that didn't pass
    % (if combinemultiple == false, this should have no effect because
    % there's only one session)
    
    if (isempty(keepers))
        
        subject_kill_list(end+1) = subject_index;
        
    else
        
        % remove these sessions
        
        temp = aap.acq_details.subjects(subject_index).seriesnumbers;
        aap.acq_details.subjects(subject_index).seriesnumbers = temp(keepers);

        temp = aap.acq_details.subjects(subject_index).mriname;
        aap.acq_details.subjects(subject_index).mriname = temp(keepers);
        
    end
    
    
end % loop over subjects

% delete any subjects on the kill list

if ~isempty(subject_kill_list)
    aap.acq_details.subjects(subject_kill_list) = [];
end

% done!

end
 


%-------------------------------------------------------------------------
function pass = eval_BIDSboldfilter(aap,value,op,test1,test2)
%-------------------------------------------------------------------------

% evaluate: value-op-test1 or (special case): test1 < value < test2

pass = false;

switch op

   case '='

       % "=" might indicate numeric equality or string equality 

       if ischar(value)
           if strcmp(value,test1); pass = true; end
       else
           if (value==test1); pass = true; end
       end

   case '>'

       if (value > test1)
           pass = true;
       end

   case '<'
       
       if (value < test1)
           pass = true;
       end
       
    case '><'
        
        if (value > test1 && value < test2)
            pass = true;
        end

    otherwise

       aas_log(aap,true,sprintf('WARNING: Unimplemented bold filter op. Halting.'));
       pass = false;

end % switch
end % eval_BIDSboldfilter


%-------------------------------------------------------------------------
function [outfname, runstr] = retrieve_file(fname)
%-------------------------------------------------------------------------

% fully specified fname with fullpath
outfname = {''};

% pre-process fname
[p, f, e] = fileparts(fname);
tags = regexp(f,'_','split');
% suffix
if isempty(regexp(f(end),'[0-9]', 'once')) % last entry is suffix
    e = [tags{end} e];
    tags(end) = [];
end
ntags = 1:numel(tags);
f = list_index(f,1,ntags);
% suffix
itags = [cell_index(tags,'acq'), cell_index(tags,'run')]; itags(~itags) = [];
runstr = '';
for t = itags
    runstr = [runstr '_' tags{t}];
end

% search
itags = [cell_index(tags,'ses'), cell_index(tags,'sub')]; itags(~itags) = [];
itagrun = cell_index(tags,'run');

if exist(fname,'file'), outfname(end+1) = {fname}; end
if itagrun
    fname = fullfile(fullfile(p,[list_index(f,1,ntags(ntags~=itagrun),0) e]));
    if exist(fname,'file'), outfname(end+1) = {fname}; end
end

p = fileparts(p);
for t = itags
    % one level up
    p = fileparts(p); ntags(t) = [];
    fname = fullfile(fullfile(p,[list_index(f,1,ntags,0) e]));
    if exist(fname,'file'), outfname(end+1) = {fname}; end
    if itagrun
        fname = fullfile(fullfile(p,[list_index(f,1,ntags(ntags~=itagrun),0) e]));
        if exist(fname,'file'), outfname(end+1) = {fname}; end
    end
end
if exist(outfname{end},'file'), outfname(1) = []; end

end
