function aap = aa_convert_subjects(aap)
aap0 = aarecipe('aap_parameters_defaults.xml','aap_tasklist_fmri.xml');
fields=fieldnames(aap0.schema.acq_details.subjects);
fields(strcmp(fields,'ATTRIBUTE')) = [];
for field=fields'
    newsubj.(field{1})={[]};
end
fields(strcmp(fields,'subjname')) = [];
newsubj.subjname = '';

newsubj(1:numel(aap.acq_details.subjects)) = newsubj;

oldsubj = aap.acq_details.subjects;
for s = 1:numel(oldsubj)
    for field=fields'
        if isfield(oldsubj(s),field{1})
            newsubj(s).(field{1})={oldsubj(s).(field{1})};
        else
            newsubj(s).(field{1})={[]};
        end
    end
end

aap.acq_details.subjects = newsubj;
aap.aap_beforeuserchanges.acq_details.subjects = newsubj;
for s = 1:numel(aap.acq_details.subjects)
    aap = aamod_evaluatesubjectnames(aap,'doit',s);
end
aap.internal.aap_initial.acq_details.subjects = aap.acq_details.subjects;