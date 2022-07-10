classdef testExamples < matlab.unittest.TestCase
    properties (TestParameter)
        dataset = testExamples.create_dataset_parameter()
    end

    methods(Test, TestTags = {'Small'})
        %% aa_downloaddemo datasets
        function downloaddemo_datasets_test(testCase)
            % Test that aa_downloaddemo datasets returns a non-empty struct
            % array, with at least a field 'ID'
            
            retval = jsondecode(fileread('datasets.json'));

            % Verification
            msg = "Empty return value from datasets.json";
            testCase.verifyNotEmpty(retval, msg);

            msg = "Values returned by datasets.json do not have a field 'ID'";
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
            switch dataset.ID
                case {'ds000114'}
                    aa_downloaddemo(aap, dataset.ID, 'sub-01')
                case {'ds002737'}
                    aa_downloaddemo(aap, dataset.ID, 'sub-01/ses-03')
                case {'LEMON_EEG'}
                    aa_downloaddemo(aap, dataset.ID, 'sub-032301')
                case {'LEMON_MRI'}
                    aa_downloaddemo(aap, dataset.ID, 'sub-032301/ses-01/anat')
                otherwise
                    aa_downloaddemo(aap, dataset.ID)
            end

            % Verification
            msg = sprintf("Directory %s does not exist after running aa_downloaddemo", aap.directory_conventions.rawdatadir);
            testCase.verifyTrue(exist(aap.directory_conventions.rawdatadir, "dir") > 0, msg)

            msg = sprintf("Directory %s is empty after running aa_downloaddemo", aap.directory_conventions.rawdatadir);
            % Directory is empty when only . and .. entries in dir listing
            testCase.verifyFalse(length(dir(aap.directory_conventions.rawdatadir)) < 3, msg)

            % Cleanup: remove temp folder
            rmdir(aap.directory_conventions.rawdatadir, 's')
        end
    end

    methods(Test, TestTags = {'Medium', 'Minimal_install'})

        %% tutorial_1
        function tutorial_1_test(testCase)
            % Test that tutorial_1 runs out-of-the-box

            tutorial_1_SPM_CH30()

            % Verification
            % TODO: How to verify correct run, apart from no-error?
        end

    end

    methods(Static)
        %% create_dataset_parameter
        function dataset_par = create_dataset_parameter()
            % Create single struct, with as fieldnames the dataset IDs.
            % That way, those field names will be the value names in the
            % test reporting
            dataset_par = struct();
            datasets = jsondecode(fileread('datasets.json'));
            for i = 1:length(datasets)
                dataset_par.(datasets(i).ID) = datasets(i);
            end
        end
    end

end
