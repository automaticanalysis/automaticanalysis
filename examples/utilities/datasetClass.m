classdef datasetClass
    properties
        ID
        URL
        type
        subset = {}
    end
    
    properties (Access=private)
        tmpdir = tempname
    end
    
    methods
        function obj = datasetClass(ID,URL,type)
            obj.ID = ID;
            obj.URL = URL;
            obj.type = type;
        end
        
        function obj = set.subset(obj,value)
            if ~iscell(value), value = strsplit(value,':'); end
            obj.subset = value;
        end
        
        function download(obj,demodir)
            % Download and unpack the data to a temp dir first
            if ~strcmp(obj.type,'AWS')
                tgz_filename = [tempname obj.type];
                options = weboptions; options.CertificateFilename = ('');
                tgz_filename = websave(tgz_filename, obj.URL, options);
            end
            switch obj.type
                case {'.zip'}
                    unzip(tgz_filename, obj.tmpdir);
                case {'.tar.gz', '.tar'}
                    untar(tgz_filename, obj.tmpdir);
                case {'AWS'}
                    if aas_shell('which aws',true,false)
                        aas_log([], true, 'AWS CLI is not installed. See <a href="https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html">AWS CLI Installation</a>' );
                    end
                    % Retrieve common files (if exist)
                    aas_shell(sprintf('aws s3 cp %s %s/%s --quiet --recursive --exclude "*" --include "task-*" --no-sign-request',obj.URL,obj.tmpdir,obj.ID));
                    aas_shell(sprintf('aws s3 cp %s %s/%s --quiet --recursive --exclude "*" --include "dwi*" --no-sign-request',obj.URL,obj.tmpdir,obj.ID));

                    if isempty(obj.subset) % obtain the whole dataset
                        aas_shell(sprintf('aws s3 cp %s %s/%s --quiet --recursive --no-sign-request',obj.URL,obj.tmpdir,obj.ID));
                    else % subset
                        for p = obj.subset
                            aas_shell(sprintf('aws s3 cp %s/%s %s/%s/%s --quiet --recursive --no-sign-request',obj.URL,p{1},obj.tmpdir,obj.ID,p{1}));
                        end
                    end
                otherwise
                    aas_log([], true, ['ERROR: Developer error, unknown dataset filetype used for downloaddemo dataset: ' dataset.filetype]);
            end
            if ~strcmp(obj.type,'AWS'), delete(tgz_filename); end
            
            if nargin > 1, obj.postprocessing(demodir); end
        end
        
        function postprocessing(obj,demodir)
            movefile(fullfile(obj.tmpdir, obj.ID, '*'), demodir);
        end
        
    end
end