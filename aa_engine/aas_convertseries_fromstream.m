% Automatic Analysis DICOM to NIFTI convertor
% Uses SPM to convert from DICOM to NIFTI-1
%  i=subject number
%  seriesnum=series number
%  outputpath=place to write NIFTI
% Rhodri Cusack MRC CBU, Cambridge 2005

% One of:
%  aas_convertseries_fromstream(aap,inputstream,[,outputpathsuffix])
%  aas_convertseries_fromstream(aap,i,inputstream,[,outputpathsuffix])
%  aas_convertseries_fromstream(aap,i,j,inputstream,[,outputpathsuffix])
%  streamname=inputstreamname
%  i=subject number; j=session number
%  seriesnum=series number
%  outputpath=place to write NIFTI
% Rhodri Cusack MRC CBU, Cambridge 2005
%
function [aap out_allechoes dicomheader subdirs]=aas_convertseries_fromstream(aap,varargin)
v=varargin;
% if (length(v)>3 && ischar(v{end-2}))
%     outputpathsuffix=v{end};
%     v(end)=[];
% else
%     outputpathsuffix='';
% end

inputstream=v{end};
v(end)=[];

% switch(length(v))
%     case 0
%         outputpath=aas_getstudypath(aap);
%     case 1
%         i=varargin{1};
%         outputpath=aas_getsubjpath(aap,i);
%     case 2
%         i=varargin{1};
%         j=varargin{2};
%         outputpath=aas_getsesspath(aap,i,j);
%     otherwise
%         aas_log(aap,true,'internal error - wrong number of arguments to aas_convertseries_fromstream.');
% end

currdir=pwd;

dicomdata=aas_getfiles_bystream(aap,v{:},inputstream);

% Make a list of all the subdirectories within this dicomdata path (written
% for multi-echoes, but probably useful otherwise) as these need to be done
% in separate chunks
subdir_index=[];
subdirs={};
for fileind=1:size(dicomdata,1)
    pth=fileparts(dicomdata(fileind,:));
    s=find(strcmp(pth,subdirs));
    if (isempty(s))
        subdirs=[subdirs pth];
        s=length(subdirs);
    else
        s=s(1);
    end
    subdir_index(fileind)=s;
end

if (isempty(dicomdata))
    aas_log(aap,1,sprintf('Did not find a dicom series called %s',dicomdirsearchpth));
end

if (length(dicomdata)==0)
    aas_log(aap,1,sprintf('Did not find any dicom data (%s) in %s',aap.directory_conventions.dicomdatafilter,dicomsearchpth));
end

out_allechoes=[];
for subdirind=1:length(subdirs)
    out=[];
    outputpath_withsuffix=fullfile(subdirs{subdirind}); %,outputpathsuffix
    aap=aas_makedir(aap,outputpath_withsuffix);
    cd(outputpath_withsuffix);
    dicomdata_subdir=dicomdata(subdir_index==subdirind,:);
    % Limit number of volumes read in at a time
    chunksize_volumes=16;
    k=1;
     
    while (k<=size(dicomdata_subdir,1))
        %    fprintf('***New chunk\n');
        oldAcquisitionNumber=-1;
        thispass_numvolumes=0;
        DICOMHEADERS=[];
        while (k<=size(dicomdata_subdir,1))
            
            tmp = spm_dicom_headers(dicomdata_subdir(k,:));
            
            try
                % [AVG] Let us use get the TR and save it to the DICOMHEADERS
                if k == 1
                    infoD = dicominfo(dicomdata_subdir(k,:));
                    
                    if strcmp(infoD.MRAcquisitionType, '3D')
                        % In 3D sequence we can find a Private field
                        % Works for Siemens scanners (not tested elsewhere)
                        str =  infoD.Private_0029_1020;
                        xstr = char(str');
                        n = findstr(xstr, 'sWiPMemBlock.adFree[8]');
                        if isempty(n)
                            error('Could not find TR in the DICOM header!')
                        end
                        [junk, r] = strtok(xstr(n:n+100), '=');
                        TR = str2double(strtok(strtok(r, '=')));
                    else
                        % If 2D sequence...
                        TR = infoD.RepetitionTime;
                    end
                    fprintf('Sequence has a TR of %.1f ms\n', TR);
                end
                
                % [AVG] Add the TR to the DICOMHEADERS explicitly before
                % saving (and in seconds!)
                tmp{1}.volumeTR = TR/1000;
            catch
            end
            
            DICOMHEADERS=[DICOMHEADERS tmp];
            
            if (DICOMHEADERS{end}.AcquisitionNumber~=oldAcquisitionNumber)
                thispass_numvolumes=thispass_numvolumes+1;
                if (thispass_numvolumes>chunksize_volumes)
                    DICOMHEADERS=DICOMHEADERS(1:(end-1));
                    break
                end
            end
            oldAcquisitionNumber=DICOMHEADERS{end}.AcquisitionNumber;
            %        fprintf('Acq %d slice %d\n',oldAcquisitionNumber,k);
            k=k+1;
        end
        
        if (~exist('echonumbers','var'))
            DICOMHEADERS_selected=DICOMHEADERS;
        else
            DICOMHEADERS_selected=[];
            for l=1:length(DICOMHEADERS);
                if(any(DICOMHEADERS{l}.EchoNumbers==echonumbers))
                    DICOMHEADERS_selected=[DICOMHEADERS_selected DICOMHEADERS(l)];
                end
            end
        end
        conv=spm_dicom_convert(DICOMHEADERS_selected,'all','flat','nii');
        out=[out(:);conv.files(:)];
    end
    out_allechoes{subdirind}=unique(out);
    if strfind(inputstream, 'dicom_structural')
        % [AVG] This is to cope with a number of strucutral images, so we
        % may have the DICOM header of each of them...
        SeriesDescription = '';
        DCMnumbers = [];
        % Loop throught the dicoms to see if the SeriesDescription changes
        for l=1:length(DICOMHEADERS);
            if ~strcmp(SeriesDescription, DICOMHEADERS{l}.SeriesDescription)
                SeriesDescription = DICOMHEADERS{l}.SeriesDescription;
                DCMnumbers = [DCMnumbers l];
            end
        end
        dicomheader={DICOMHEADERS{DCMnumbers}};
    else
        dicomheader{subdirind}=DICOMHEADERS{1};
    end
end

% Single echo, no echo dimension to cell array
if (length(out_allechoes)==1)
    out_allechoes=out_allechoes{1};
end

cd (currdir);