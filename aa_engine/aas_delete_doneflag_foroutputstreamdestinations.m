% Recursively delete the done flags of all stages that take output from
% stage k
% aap=aas_delete_doneflag_foroutputstreamdestinations(aap,k[,i[,j]])

function aap=aas_delete_doneflag_foroutputstreamdestinations(aap,k,varargin)

for destind=1:length(aap.internal.outputstreamdestinations{k})
    deststage=aap.internal.outputstreamdestinations{k}(destind);
    aap=aas_delete_doneflag(aap,deststage,varargin{:});
    aap=aas_delete_doneflag_foroutputstreamdestinations(aap,deststage,varargin{:});
end;