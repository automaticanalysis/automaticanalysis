classdef (TestTags = {'Large'}) ...
        testUseCases < matlab.unittest.TestCase
    %TESTUSECASES Tests that represent a full use case.

    properties (Constant)
        testdir = fileparts(mfilename("fullpath"));
    end

    properties (TestParameter)
        testname = testUseCases.get_testscripts()
    end

    properties
        inputargs
    end

    methods (TestClassSetup)
        function classSetup(testCase)
            testCase.inputargs = testUseCases.pass_inputargs('get');
        end
    end

    methods (Test)
        function use_case_test(testCase, testname)
            savedir = pwd;
            cd(testCase.testdir);

            testUseCases.runit(testname, testCase.inputargs.parameterfile, testCase.inputargs.deleteprevious, testCase.inputargs.haltonerror, testCase.inputargs.wheretoprocess);

            cd(savedir);
        end
    end

    methods (Static)
        function testnames = get_testscripts()
            tests = dir(fullfile(testUseCases.testdir,'aatest_*.m'));
            testnames = {tests.name};
        end

        function runit(scriptname, parameterfile, deleteprevious, haltonerror, wheretoprocess)

            func = str2func(strrep(scriptname,'.m',''));

            if  haltonerror

                % run script normally;
                % function halts on aa error

                func(parameterfile, deleteprevious, wheretoprocess);

                % if encapsulated returns, this script passed

                fprintf('\npass - %s\n',scriptname);

            else

                % run script inside a try-catch;
                % function flags any error and returns
                % optionally dropping into qsub_debug

                try

                    func(parameterfile, deleteprevious, wheretoprocess);
                    fprintf('\npass - %s\n',scriptname);

                catch err

                    fprintf('FAIL - %s\n',scriptname);

                    if strcmp(wheretoprocess,'localsingle')
                        testUseCases.reporterror(scriptname,err);
                    else

                        try
                            aaq_qsub_debug;
                        catch err
                            testUseCases.reporterror(scriptname,err);
                        end

                        % catch cases when aaq_qsub_debug quietly
                        % returns (e.g., failed to submitjob)

                        keyboard;

                    end

                end

            end

        end

        function reporterror(scriptname,err)
            msg = sprintf('%s had an error: %s\n',scriptname,err.message);
            for e = 1:numel(err.stack)
                % Stop tracking to internal
                if strfind(err.stack(e).file,'distcomp'), break, end
                msg = [msg sprintf('<a href="matlab: opentoline(''%s'',%d)">in %s (line %d)</a>\n', ...
                    err.stack(e).file, err.stack(e).line,...
                    err.stack(e).file, err.stack(e).line)];
            end
            fprintf(msg)
        end

        function output = pass_inputargs(action, varargin)
            persistent inputargs
            output = true;
            switch action
                case 'get'
                    output = inputargs;
                case 'set'
                    inputargs = varargin{1};
            end
        end

    end
end

