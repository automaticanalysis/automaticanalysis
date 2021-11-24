classdef testExamples < matlab.unittest.TestCase
    properties (TestParameter)
        dataset = num2cell(aa_downloaddemo_datasets())
    end

    methods(Test, TestTags = {'Small'})
        %% aa_downloaddemo_datasets
        function downloaddemo_datasets_test(testCase)
            % Test that aa_downloaddemo_datasets returns a non-empty struct
            % array, with at least a field 'ID'

            retval = aa_downloaddemo_datasets();

            % Verification
            msg = "Empty return value from aa_downloaddemo_datasets";
            testCase.verifyNotEmpty(retval, msg);

            msg = "Values returned by aa_downloaddemo_datasets do not have a field 'ID'";
            testCase.verifyTrue(isfield(retval, 'ID'), msg)
        end

    end

    methods(Test, TestTags = {'Medium'})

        %% aa_downloaddemo
        % Note: although the contents of this test are 'Small', it is tagged
        % as 'Medium' because the actual download can take a while.
        function downloaddemo_test(testCase, dataset)
            % Test that all predefined data downloads work correctly

            % Use tempname to ensure directory does not exist yet, so can
            % verify at the end that exists and non-empty
            aap.directory_conventions.rawdatadir = tempname;

            % Add minimal fields to use aap in aas_log
            % TODO: have a minimal_aap constructor somewhere, at least for
            % use in testing?
            aap.options.verbose = 2;
            aap.options.email = '';
            aap.gui_controls.usecolouroutput = false;

            % Run the aa_downloaddemo
            aa_downloaddemo(aap, dataset.ID)

            % Verification
            msg = sprintf("Directory %s does not exist after running aa_downloaddemo", aap.directory_conventions.rawdatadir);
            testCase.verifyTrue(exist(aap.directory_conventions.rawdatadir, "dir") > 0, msg)

            msg = sprintf("Directory %s is empty after running aa_downloaddemo", aap.directory_conventions.rawdatadir);
            % Directory is empty when only . and .. entries in dir listing
            testCase.verifyFalse(length(dir(aap.directory_conventions.rawdatadir)) < 3, msg)

            % Cleanup: remove temp folder
            rmdir(aap.directory_conventions.rawdatadir, 's')
        end

        %% tutorial_1_aa_setup
        function tutorial_1_test(testCase)
            % Test that tutorial_1_aa_setup runs out-of-the-box

            tutorial_1_aa_setup()

            % Verification
            % TODO: How to verify correct run, apart from no-error?
        end
    end
end