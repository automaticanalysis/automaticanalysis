function aa_report_export(studyPath, target)

    mediaDir = fullfile(target,'media');

    aas_makedir([],target);
    aas_makedir([],mediaDir);

    load(fullfile(studyPath,'aap_parameters_reported.mat'),'aap');
    oldRoot = fullfile(aap.acq_details.root,aap.directory_conventions.analysisid);

    copyfile(fullfile(oldRoot,'aa_*'),target);
    copyfile(fullfile(oldRoot,'aap_*'),target);
    for fn = reshape(aap.report.dependency,1,[])
        copyfile(fn{1},mediaDir);
    end

    reportFields = fieldnames(aap.report);

    for repDir = reshape(reportFields(contains(fieldnames(aap.report),'dir')),1,[])
        aas_makedir([],fullfile(strrep(aap.report.(repDir{1}),oldRoot,target)));
    end

    % top-level HTMLs
    topHMLSs = {'html_main' 'html_S00' 'html_moco' 'html_reg' 'html_C00'};

    % HTMLs in subfolders
    subHTMLs = reportFields(cellfun(@(f) ~contains(f,{'S00' 'C00'}) & ~isempty(regexp(f,'(S[0-9]{2})|(C[0-9]{2})', 'once')), reportFields))';

    for html = [topHMLSs subHTMLs]
        switch html{1}
            case topHMLSs
                relTarget = '.';
                relMedia = './media/';
            case subHTMLs
                relTarget = '..';
                relMedia = '../media/';
        end

        content = strrep(fileread(aap.report.(html{1}).fname),'\','/');
        content = regexprep(content,['(?<=href=")' strrep(oldRoot,'\','/')],relTarget);
        content = regexprep(content,'(?<=src=")[a-zA-Z0-9-_:\\/]*(?=diagnostic)',relMedia);

        newFn = strrep(aap.report.(html{1}).fname,oldRoot,target);
        fid = fopen(newFn,'w');
        if fid == -1, aas_log(aap,true,sprintf('Failed to open %s',newFn)); end
        fprintf(fid,'%s',content);
        fclose(fid);
    end
end
